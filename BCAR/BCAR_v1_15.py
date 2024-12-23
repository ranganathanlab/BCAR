# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:10:30 2024

@author: bryan

Last debugging on 12/22/2024
"""

import gzip
import multiprocessing as mp
import sys
import argparse
import time
import math
from .needleman_wunsch import needleman_wunsch
from .build_consensus import build_consensus

#initialize generic things
twobit_basespace = {'A':0, 'C':1, 'G':2, 'T':3}
alignment_basespace = {"-": 0,'A': 1,'C': 2,'G': 3,'T': 4,'a': 5,'c': 6,'g': 7,'t': 8,'N': 9}
rev_bases =           {0: '-',1: 'A',2: 'C',3: 'G',4: 'T',5: 'a',6: 'c',7: 'g',8: 't',9: 'N'}
start_time = time.time()

alignment_time = time.time()

def open_fastq_file(filename):
    """Opens a FASTQ file, handling gzip compression if necessary."""
    if filename.endswith('.fastq.gz'):
        return gzip.open(filename, 'rt')  # Open gzipped FASTQ in text mode
    elif filename.endswith('.fastq'):
        return open(filename, 'r')  # Open regular FASTQ file
    else:
        raise ValueError("Unsupported file format. File must end with .fastq or .fastq.gz")

def parse_arguments():
    """Parses command-line arguments using argparse."""
    parser = argparse.ArgumentParser(description="Process FASTQ files and generate a barcode-sequence map.")

    parser.add_argument('--fwd', dest='fwd_fastqs', type=str, nargs='+', required=True,
                        help="List of .fastq or .fastq.gz files containing reads")
    parser.add_argument('--rev', dest='rev_fastqs', type=str, nargs='+', required=True,
                        help="List of .fastq or .fastq.gz files containing reads")
    parser.add_argument('--out1', dest='output_fastq_fwd', type=str, default="./out.fwd.fastq",
                        help=".fastq file that will contain the consensus forward reads")
    parser.add_argument('--out2', dest='output_fastq_rev', type=str, default="./out.rev.fastq",
                        help=".fastq file that will contain the consensus reverse reads")
    parser.add_argument('--BC-start', dest='bc_start', type=int, default=0,
                        help="Position in the sequence where the barcode begins")
    parser.add_argument('--BC-len', dest='bc_len', type=int, default=15,
                        help="Length of the barcode")
    parser.add_argument('-t', '--threads', dest='threads', type=int, default=0,
                        help="Number of threads to use (if left empty, will attempt to use max number available.)")
    parser.add_argument('--min-qscore', dest='min_qscore', type=int, default=20,
                        help="Reads are trashed if a base in the barcode or mapped region is below this quality score")
    parser.add_argument('--min-count', dest='min_count', type=int, default=1,
                        help="Minimum number of times to observe a barcode before it is reported")
    parser.add_argument('-a', '--align', action='store_true',
                        help="Whether to align reads prior to calling consensus (useful for indel-prone sequencers like NanoPore / Aviti)")
    return parser.parse_args()

def seq_to_twobit(seq, encode_size = 4):
    # converts a string of bases into an integer for data compression
    # only accepts A,C,G,T, so filter reads before getting here
    int_list = [twobit_basespace[seq[i]]*encode_size**i for i in range(len(seq))]
    #add capping "C" to prevent trailing A's from getting ignored (since they would add zero)
    int_list.append(encode_size**len(seq))
    return sum(int_list)

def twobit_to_seq(seq_int, encode_size = 4):
    # converst a two-bit integer representing a sequence back into a string of bases
    seq = ''.join([rev_bases[b] for b in twobit_to_list(seq_int, encode_size=encode_size)])
    return seq

def twobit_to_list(seq_int, encode_size=4): #encoded in ACTG -> 0123
    # converts two-bit encoded integer to base-10 encoded list
    seq_len = max([int(math.log(seq_int, encode_size)),0])
    seq_list = []
    for i in range(seq_len, -1, -1): #backwards through sequence
        seq_list.append(seq_int//encode_size**i+1) # +1 to encode to ACTG -> 1234
        seq_int = seq_int % encode_size**i
    seq_list = seq_list[::-1] #reverse list
    return seq_list[:-1] #remove capping "C"

def multi_align_lists(ali_lists):
    # determine max possible sequence length after inserting gaps
    seq_len = max([len(seq) for seq in ali_lists])
    total_inserts = sum([sum([seq.count(i) for i in [5,6,7,8]]) for seq in ali_lists])
    seq_len += total_inserts
    
    N = len(ali_lists) #number of sequences

    # insert gaps wherever insertions relative to the reference appear in other seqs
    for i in range(seq_len):
        bases = [0 for _ in range(N)]
        for j in range(N):
            if i < len(ali_lists[j]):
                bases[j] = ali_lists[j][i]
        if max(bases) <= 4: #no insertions
            continue
        else:
            for j in range(N):
                if bases[j] <= 4:
                    ali_lists[j].insert(i, 0)

    return ali_lists

def make_barcode_printables(list_of_twobit_bcs, fwd_list_of_lists, rev_list_of_lists, min_count, align_flag):
    conversion_time = 0.
    alignment_time  = 0.
    consensus_time  = 0.
    
    
    fwd_results_list = [] #where strings to print will go
    rev_results_list = []
    for i in range(len(list_of_twobit_bcs)):

        # check count and skip if below minimum
        count = len(fwd_list_of_lists[i])
        if count < min_count:
            continue
        
        #if only one read, no need to align or build consensus
        if count == 1:
            bc_seq = twobit_to_seq(list_of_twobit_bcs[i])
            fwd_seq = twobit_to_seq(fwd_list_of_lists[i][0])
            rev_seq = twobit_to_seq(rev_list_of_lists[i][0])
            fwd_results_list.append("@fwd_consensus_bc_read;bc=%s;count=1\n%s\n+\n%s\n" % (
                bc_seq,
                fwd_seq,
                ''.join(["4" for q in range(len(fwd_seq))]))
                )
            rev_results_list.append("@rev_consensus_bc_read;bc=%s;count=1\n%s\n+\n%s\n" % (
                bc_seq,
                rev_seq,
                ''.join(["4" for q in range(len(rev_seq))]))
                )
            continue
        
        #else, do full processing
        # convert out of integer space
        temp_time = time.time()
        bc_seq = twobit_to_seq(list_of_twobit_bcs[i])
        fwd_seqs = [twobit_to_list(seq) for seq in fwd_list_of_lists[i]]
        rev_seqs = [twobit_to_list(seq) for seq in rev_list_of_lists[i]]
        conversion_time += (time.time()-temp_time)
        # align if requested
        
        if align_flag:
            temp_time = time.time()
            fwd_ali_seqs = [needleman_wunsch(fwd_seqs[0], seq) for seq in fwd_seqs]
            fwd_ali_seqs = multi_align_lists(fwd_ali_seqs)
            rev_ali_seqs = [needleman_wunsch(rev_seqs[0], seq) for seq in rev_seqs]
            rev_ali_seqs = multi_align_lists(rev_ali_seqs)
            alignment_time += (time.time()-temp_time)
        else:
            fwd_ali_seqs = fwd_seqs
            rev_ali_seqs = rev_seqs
        
        # make consensus sequences
        temp_time = time.time()
        fwd_cons_seq, fwd_cons_q = build_consensus(fwd_ali_seqs)
        rev_cons_seq, rev_cons_q = build_consensus(rev_ali_seqs)
        consensus_time += (time.time()-temp_time)
        
        # format to fastq strings
        fwd_results_list.append("@fwd_consensus_bc_read;bc=%s;count=%s\n%s\n+\n%s\n" % (
            bc_seq,
            count,
            ''.join([rev_bases[b] for b in fwd_cons_seq]),
            ''.join([chr(q + 33) for q in fwd_cons_q]))
            )
        rev_results_list.append("@rev_consensus_bc_read;bc=%s;count=%s\n%s\n+\n%s\n" % (
            bc_seq,
            count,
            ''.join([rev_bases[b] for b in rev_cons_seq]),
            ''.join([chr(q + 33) for q in rev_cons_q]))
            )
    #sys.stderr.write("spent %s seconds on conversion, %s seconds on alignment, and %s seconds on consensus building\n" % (conversion_time, alignment_time, consensus_time))
    return fwd_results_list, rev_results_list

def read_fastq_chunk(open_fastq_file, chunk_size=50000):
    # read chunk_size reads into memory from a fastq. Does not trim newlines.
    # also checks for end of file
    eof = False
    reads = []
    qscores = []
    for i in range(chunk_size):
        header = open_fastq_file.readline()
        if not header: #check for eof
            eof = True
            break
        reads.append(open_fastq_file.readline())
        _ = open_fastq_file.readline()
        qscores.append(open_fastq_file.readline())
    return reads, qscores, eof

def process_read(fwd_read, fwd_qscore, rev_read, rev_qscore, min_qscore, bc_start, bc_len, bc_end):
    # Add quality filters here
    # check qscores
    for q in fwd_qscore[:-1]:
        if ord(q) < min_qscore:
            return None, None, None
    for q in rev_qscore[:-1]:
        if ord(q) < min_qscore:
            return None, None, None
    # Make sure you're not trying to grab a barcode from a read that's too short
    if len(fwd_read) <= bc_len or len(rev_read) <= bc_len:
        return None, None, None
    
    bc_int  = seq_to_twobit(fwd_read[bc_start:bc_end])
    fwd_int = seq_to_twobit(fwd_read[:-1])
    rev_int = seq_to_twobit(rev_read[:-1])
    return bc_int, fwd_int, rev_int

def process_minibatch(fwd_reads, fwd_qscores, rev_reads, rev_qscores, min_qscore, bc_start, bc_len, bc_end):
    # generally faster with mulitprocessing to hand each thread many reads at a time
    minibatch_fwd_map = {}
    minibatch_rev_map = {}
    for i in range(len(fwd_reads)):
        bc_int, fwd_int, rev_int = process_read(fwd_reads[i], fwd_qscores[i], rev_reads[i], rev_qscores[i], min_qscore, bc_start, bc_len, bc_end)
        if bc_int:
            if minibatch_fwd_map.get(bc_int) == None:
                minibatch_fwd_map[bc_int] = [fwd_int]
                minibatch_rev_map[bc_int] = [rev_int]
            else:
                minibatch_fwd_map[bc_int].append(fwd_int)
                minibatch_rev_map[bc_int].append(rev_int)
    return minibatch_fwd_map, minibatch_rev_map

def main():
    args = parse_arguments()
    # initialize dataset-specific things
    mb_size = 1000                                   #size of minibatch handed to each process
    chunk_size = 1000000                             #size of chunk to read from file (in reads)
    args.min_qscore +=33                      
    bc_end  = args.bc_start + args.bc_len
    if args.threads == 0:
        args.threads = mp.cpu_count()
    
    # Read over fastq file and call variants
    bc_fwd_map = {} #will contain all barcodes and occupy most used memory
    bc_rev_map = {} #matched map for read 2
    with mp.Pool(processes=args.threads) as pool:
        for file_num in range(len(args.fwd_fastqs)):
            with open_fastq_file(args.fwd_fastqs[file_num]) as open_fastq_fwd, open_fastq_file(args.rev_fastqs[file_num]) as open_fastq_rev:
                eof_check = False
                chunk_count = 0

                while not eof_check:
                    # Read a chunk into memory
                    fwd_read_chunk, fwd_qscore_chunk, eof_check = read_fastq_chunk(open_fastq_fwd, chunk_size=chunk_size)
                    rev_read_chunk, rev_qscore_chunk, eof_check = read_fastq_chunk(open_fastq_rev, chunk_size=chunk_size)

                    # Break chunk into minibatches and send to processes
                    results = pool.starmap(
                        process_minibatch,
                        [
                            (fwd_read_chunk[i:i + mb_size], fwd_qscore_chunk[i:i + mb_size], 
                             rev_read_chunk[i:i + mb_size], rev_qscore_chunk[i:i + mb_size], 
                             args.min_qscore, args.bc_start, args.bc_len, bc_end)
                            for i in range(0, len(fwd_read_chunk), mb_size)
                        ]
                    )

                    # Add sequences to main bc_map
                    for mb_fwd_map, mb_rev_map in results:
                        for bc in mb_fwd_map:
                            if bc_fwd_map.get(bc) is None:
                                bc_fwd_map[bc] = mb_fwd_map[bc]
                                bc_rev_map[bc] = mb_rev_map[bc]
                            else:
                                bc_fwd_map[bc] += mb_fwd_map[bc]
                                bc_rev_map[bc] += mb_rev_map[bc]

                    # Keep tabs on progress
                    chunk_count += 1
                    if not eof_check:
                        sys.stderr.write(f"Processed {chunk_count * chunk_size} reads from {args.fwd_fastqs[file_num]}\n")

    sys.stderr.write("Called %s barcodes in %s seconds\n" % (len(bc_fwd_map), round(time.time()-start_time,2)))
    sys.stderr.write("Writing to file...\n")
    
    # Convert dictionaries to matched-order lists
    bc_list  = list(bc_fwd_map.keys())
    fwd_list = []
    rev_list = []
    for i in range(len(bc_list)):
        fwd_list.append(bc_fwd_map.pop(bc_list[i]))
        rev_list.append(bc_rev_map.pop(bc_list[i]))

    # Print to file
    with open(args.output_fastq_fwd,'w+') as fwd_out, open(args.output_fastq_rev,'w+') as rev_out:
        with mp.Pool(processes=args.threads) as pool:
            batch_size = 100
            results = []
            # Set up solve
            for i in range(0, len(bc_list), batch_size):
                        args_tuple = (bc_list[i:i + batch_size],
                                      fwd_list[i:i + batch_size],
                                      rev_list[i:i + batch_size],
                                      args.min_count,
                                      args.align)
                        results.append(pool.apply_async(make_barcode_printables, args_tuple))
            for result in results:
                res = result.get()  # Collect results
                if res is not None:
                    for j in range(batch_size):
                        try:
                            fwd_out.write(res[0][j])
                            rev_out.write(res[1][j])
                        except IndexError:
                            break
    sys.stderr.write("finished in %s seconds\n" % round(time.time()-start_time,2))

if __name__ == "__main__":
    main()

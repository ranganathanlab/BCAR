# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:10:30 2024

@author: bryan
"""

import gzip
import multiprocessing as mp
import sys
import argparse
import time
import math
from needleman_wunsch import needleman_wunsch

#initialize generic things
twobit_basespace = {'A':0, 'C':1, 'G':2, 'T':3}
alignment_basespace = {"-": 0,'A': 1,'C': 2,'G': 3,'T': 4,'a': 5,'c': 6,'g': 7,'t': 8,'N': 9}
rev_bases =           {0: '-',1: 'A',2: 'C',3: 'G',4: 'T',5: 'a',6: 'c',7: 'g',8: 't',9: 'N'}
start_time = time.time()

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
    parser.add_argument('--BC-start', dest='bc_start', type=int, default=196,
                        help="Position in the sequence where the barcode begins")
    parser.add_argument('--BC-len', dest='bc_len', type=int, default=15,
                        help="Length of the barcode")
    parser.add_argument('--min-qscore', dest='min_qscore', type=int, default=20,
                        help="Reads are trashed if a base in the barcode or mapped region is below this quality score")
    parser.add_argument('--min-agree', dest='min_agreement', type=float, default=0.7,
                        help="Fraction of bases at each position in the consensus that must agree, else the barcode is trashed")
    parser.add_argument('--min-count', dest='min_count', type=int, default=1,
                        help="Minimum number of times you see a barcode before you map it")

    return parser.parse_args()

def seq_to_twobit(seq, encode_size = 4):
    #converts a string of bases into an integer for data compression
    #only accepts A,C,G,T, so filter reads before getting here
    int_list = [twobit_basespace[seq[i]]*encode_size**i for i in range(len(seq))]
    return sum(int_list)

def twobit_to_seq(seq_int, encode_size = 4):
    #converst a two-bit integer representing a sequence back into a string of bases
    seq = ''.join([rev_bases[b] for b in twobit_to_list(seq_int, encode_size=encode_size)])
    return seq

def twobit_to_list(seq_int, encode_size=4): #encoded in ACTG -> 0123
    #converts two-bit encoded integer to base-10 encoded list
    seq_len = int(math.log(seq_int, encode_size))
    seq_list = []
    for i in range(seq_len, -1, -1): #backwards through sequence
        seq_list.append(seq_int//encode_size**i+1)
        seq_int = seq_int % encode_size**i
    return seq_list[::-1] #encoded in ACTG -> 1234

def multi_align_lists(ali_lists):
    #determine max possible sequence length after inserting gaps
    seq_len = max([len(seq) for seq in ali_lists])
    total_inserts = sum([sum([seq.count(i) for i in [5,6,7,8]]) for seq in ali_lists])
    seq_len += total_inserts
    
    N = len(ali_lists) #number of sequences

    #insert gaps wherever insertions relative to the reference appear in other seqs
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

    #return [''.join([rev_bases[b] for b in ali_lists[j]]) for j in range(N)]
    return ali_lists

def create_consensus(list_of_seq_ints, min_agree = 0.7, min_count = 1):
    if len(list_of_seq_ints) < min_count:
        sys.stderr.write("too few sequences for a barcode\n")
        return None
    #align each seq to reference
    seqs = [twobit_to_list(seq) for seq in list_of_seq_ints]
    ali_seqs = [needleman_wunsch(seqs[0], seq) for seq in seqs]
    ali_seqs = multi_align_lists(ali_seqs)
    
    #construct consensus sequence
    max_len = max([len(seq) for seq in ali_seqs])
    consensus = []
    for i in range(max_len):
        base_counts = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0] #base-10 encoding
        for seq in ali_seqs:
            try:
                base_counts[seq[i]] += 1
            except IndexError:
                continue
        if max(base_counts) / sum(base_counts) >= min_agree:
            if base_counts[0] != max(base_counts): #don't append gaps
                consensus.append(base_counts.index(max(base_counts)))
        else:
            consensus.append(9)
    consensus = [rev_bases[b] for b in consensus]
    return ''.join(consensus)

def make_scores_for_read(consensus_read):
    qscores = ["a" for base in consensus_read]
    for i in range(len(qscores)):
        if consensus_read[i] == "N":
            qscores[i] = "J"
    return ''.join(qscores)

def make_barcode_printables(list_of_twobit_bcs, fwd_list_of_lists, rev_list_of_lists, min_agree, min_count):
    fwd_results_list = []
    rev_results_list = []
    for i in range(len(list_of_twobit_bcs)):
        bc_seq = twobit_to_seq(list_of_twobit_bcs[i])
        count = len(fwd_list_of_lists[i])
        consensus_fwd = create_consensus(fwd_list_of_lists[i], min_agree = min_agree, min_count = min_count)
        consensus_rev = create_consensus(rev_list_of_lists[i], min_agree = min_agree, min_count = min_count)
        if consensus_fwd and consensus_rev:
            fwd_results_list.append("@fwd_consensus_bc_read;bc=%s;count=%s\n%s\n+\n%s\n" % (bc_seq, count, consensus_fwd, make_scores_for_read(consensus_fwd)))
            rev_results_list.append("@rev_consensus_bc_read;bc=%s;count=%s\n%s\n+\n%s\n" % (bc_seq, count, consensus_rev, make_scores_for_read(consensus_rev)))
    return fwd_results_list, rev_results_list

def read_fastq_chunk(open_fastq_file, chunk_size=50000):
    #read chunk_size reads into memory from a fastq. Does not trim newlines.
    #also checks for end of file
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

def preprocess_read(fwd_read, fwd_qscore, rev_read, rev_qscore, min_qscore, bc_start, bc_end):
    #check qscores
    for q in fwd_qscore[:-1]:
        if ord(q) < min_qscore:
            return None, None, None
    for q in rev_qscore[:-1]:
        if ord(q) < min_qscore:
            return None, None, None
    #filter on length
    if len(fwd_read) < 50 or len(rev_read) < 50:
        return None, None, None
    #other filters can be put in here like length, etc.
    #return fwd and reverse reads with barcode
    return fwd_read[:-1], rev_read[:-1], fwd_read[bc_start:bc_end] #trim newlines

def process_read(fwd_read, fwd_qscore, rev_read, rev_qscore, min_qscore, bc_start, bc_len, bc_end):
    fwd_read, rev_read, bc = preprocess_read(fwd_read, fwd_qscore, rev_read, rev_qscore, min_qscore, bc_start, bc_end)
    if not bc:
        return None, None, None
    bc_int  = seq_to_twobit(bc)
    fwd_int = seq_to_twobit(fwd_read)
    rev_int = seq_to_twobit(rev_read)
    return bc_int, fwd_int, rev_int

def process_minibatch(fwd_reads, fwd_qscores, rev_reads, rev_qscores, min_qscore, bc_start, bc_len, bc_end):
    #generally faster with mulitprocessing to hand each thread many reads at a time
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
    #initialize dataset-specific things
    mb_size = 1000                                   #size of minibatch handed to each process
    chunk_size = 1000000                             #size of chunk to read from file (in reads)
    args.min_qscore +=33                      
    bc_end  = args.bc_start + args.bc_len
    
    # Read over fastq file and call variants
    bc_fwd_map = {} #will contain all barcodes and occupy most used memory
    bc_rev_map = {} #matched map for read 2
    with mp.Pool(processes=mp.cpu_count()) as pool:
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
        # Write to file
        num_minibatches = chunk_size // mb_size // 10
        
    sys.stderr.write("Called %s barcodes in %s seconds\n" % (len(bc_fwd_map), round(time.time()-start_time,2)))
    sys.stderr.write("Writing to file...\n")
    
    # Write to file
    batch_size = 1000 #barcodes per process to find consensus and print
    cpu_count = mp.cpu_count()
    with open(args.output_fastq_fwd,'w+') as fwd_out:
        with open(args.output_fastq_rev,'w+') as rev_out:
            with mp.Pool(processes=cpu_count) as pool:
                i=0
                bc_batches  = []
                fwd_batches = []
                rev_batches = []
                for bc in bc_fwd_map:
                    i+=1
                    bc_batches.append(bc)
                    fwd_batches.append(bc_fwd_map[bc])
                    rev_batches.append(bc_rev_map[bc])
                    if i >= (batch_size * cpu_count):
                        results= pool.starmap(
                            make_barcode_printables,
                            [(bc_batches[j:j+batch_size], fwd_batches[j:j+batch_size], rev_batches[j:j+batch_size],
                             args.min_agreement, args.min_count)
                            for j in range(0,i,batch_size)]
                        )
                        for res in results:
                            for j in range(batch_size):
                                fwd_out.write(res[0][j])
                                rev_out.write(res[1][j])
                        i=0
                        bc_batches  = []
                        fwd_batches = []
                        rev_batches = []
                #last loop around, pick up any remaining after batch size
                batch_size = len(bc_batches)//cpu_count +1
                results= pool.starmap(
                    make_barcode_printables,
                    [(bc_batches[j:j+batch_size], fwd_batches[j:j+batch_size], rev_batches[j:j+batch_size],
                     args.min_agreement, args.min_count)
                    for j in range(0,i,batch_size)]
                )
                for res in results:
                    for j in range(batch_size):
                        try:
                            fwd_out.write(res[0][j])
                            rev_out.write(res[1][j])
                        except IndexError:
                            break
                        
    sys.stderr.write("finished in %s seconds\n" % round(time.time()-start_time,2))

if __name__ == "__main__":
    main()

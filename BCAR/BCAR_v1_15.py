#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 12:45:25 2024

@author: bryan
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
twobit_basespace    = {'A':0, 'C':1, 'G':2, 'T':3}
alignment_basespace = {"-":0, 'A':1, 'C':2, 'G':3, 'T':4, 'a':5, 'c':6, 'g':7 ,'t':8, 'N':9}
rev_bases           = {0:'-', 1:'A', 2:'C', 3:'G', 4:'T', 5:'a', 6:'c', 7:'g', 8:'t', 9:'N'}
start_time = time.time()

def parse_arguments():
    """Parses command-line arguments using argparse."""
    parser = argparse.ArgumentParser(description="Process FASTQ files and generate a barcode-sequence map.")

    parser.add_argument('--fwd', dest='fwd_fastqs', type=str, nargs='+', required=True,
                        help="List of .fastq or .fastq.gz files containing reads")
    parser.add_argument('--rev', dest='rev_fastqs', type=str, nargs='+', required=False,
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

def main():
    args = parse_arguments()
    # initialize dataset-specific things
    args.min_qscore +=33                      
    bc_end  = args.bc_start + args.bc_len
    if args.threads == 0:
        args.threads = mp.cpu_count()
    if args.rev_fastqs:
        paired_mode = True
    else:
        paired_mode = False
    
    # collect barcode-read pairs into dictionary
    if paired_mode:
        bc_fwd_map, bc_rev_map = collect_barcodes_paired(
            args.fwd_fastqs, args.rev_fastqs, args.min_qscore, args.bc_start, args.bc_len, bc_end, args.threads
            )
    else:
        bc_fwd_map = collect_barcodes_single(
            args.fwd_fastqs, args.min_qscore, args.bc_start, args.bc_len, bc_end, args.threads
            )
    sys.stderr.write("Collected %s barcodes in %s seconds\n" % (len(bc_fwd_map), round(time.time()-start_time,2)))

    # Convert barcode-sequence map to ordered lists (frees memory and makes zippable)
    bc_list   = list(bc_fwd_map.keys())
    fwd_list  = [bc_fwd_map.pop(bc) for bc in bc_list]
    if paired_mode:
        rev_list  = [bc_rev_map.pop(bc) for bc in bc_list]
    
    # Simultaneously collapse barcodes and write them to file
    sys.stderr.write("Writing to file...\n")
    if paired_mode:
        write_to_file_paired(bc_list, fwd_list, rev_list, args.output_fastq_fwd, args.output_fastq_rev, args.min_count, args.align, args.threads)
    else:
        write_to_file_single(bc_list, fwd_list, args.output_fastq_fwd, args.min_count, args.align, args.threads)

    sys.stderr.write("finished in %s seconds\n" % round(time.time()-start_time,2))

#############################################
#   Called in main() / High-level fuctions  #
#############################################

def collect_barcodes_paired(fwd_fastqs, rev_fastqs, min_qscore, bc_start, bc_len, bc_end, threads):
    chunk_size = 1000000  #size of chunk to read from file (in reads)
    mb_size = 1000        #size of minibatch, aught to divide into chunk_size
    
    bc_fwd_map = {} #will contain all barcodes and occupy most used memory
    bc_rev_map = {} #matched map for read 2
    with mp.Pool(processes=threads) as pool:
        for file_num in range(len(fwd_fastqs)):
            with open_fastq_file(fwd_fastqs[file_num]) as open_fastq_fwd, open_fastq_file(rev_fastqs[file_num]) as open_fastq_rev:
                eof_check = False
                chunk_count = 0

                while not eof_check:
                    # Read a chunk into memory
                    fwd_read_chunk, fwd_qscore_chunk, eof_check = read_fastq_chunk(open_fastq_fwd, chunk_size=chunk_size)
                    rev_read_chunk, rev_qscore_chunk, eof_check = read_fastq_chunk(open_fastq_rev, chunk_size=chunk_size)

                    # Break chunk into minibatches and send to processes
                    results = pool.starmap(
                        process_minibatch_paired,
                        [
                            (fwd_read_chunk[i:i + mb_size], fwd_qscore_chunk[i:i + mb_size], 
                             rev_read_chunk[i:i + mb_size], rev_qscore_chunk[i:i + mb_size], 
                             min_qscore, bc_start, bc_len, bc_end)
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
                        sys.stderr.write(f"Processed {chunk_count * chunk_size} reads from {fwd_fastqs[file_num]}\n")
    return bc_fwd_map, bc_rev_map

def collect_barcodes_single(fwd_fastqs, min_qscore, bc_start, bc_len, bc_end, threads):
    chunk_size = 1000000  #size of chunk to read from file (in reads)
    mb_size = 1000        #size of minibatch, aught to divide into chunk_size
    
    bc_fwd_map = {} #will contain all barcodes and occupy most used memory
    with mp.Pool(processes=threads) as pool:
        for file_num in range(len(fwd_fastqs)):
            with open_fastq_file(fwd_fastqs[file_num]) as open_fastq_fwd:
                eof_check = False
                chunk_count = 0

                while not eof_check:
                    # Read a chunk into memory
                    fwd_read_chunk, fwd_qscore_chunk, eof_check = read_fastq_chunk(open_fastq_fwd, chunk_size=chunk_size)

                    # Break chunk into minibatches and send to processes
                    results = pool.starmap(
                        process_minibatch_single,
                        [
                            (fwd_read_chunk[i:i + mb_size], fwd_qscore_chunk[i:i + mb_size],
                             min_qscore, bc_start, bc_len, bc_end)
                            for i in range(0, len(fwd_read_chunk), mb_size)
                        ]
                    )

                    # Add sequences to main bc_map
                    for mb_fwd_map in results:
                        for bc in mb_fwd_map:
                            if bc_fwd_map.get(bc) is None:
                                bc_fwd_map[bc] = mb_fwd_map[bc]
                            else:
                                bc_fwd_map[bc] += mb_fwd_map[bc]

                    # Keep tabs on progress
                    chunk_count += 1
                    if not eof_check:
                        sys.stderr.write(f"Processed {chunk_count * chunk_size} reads from {fwd_fastqs[file_num]}\n")
    return bc_fwd_map

def write_to_file_paired(bc_list, fwd_list, rev_list, output_fastq_fwd, output_fastq_rev, min_count, align_flag, threads):
    min_list  = [min_count] * len(bc_list)
    flag_list = [align_flag] * len(bc_list)
    with open(output_fastq_fwd, 'w+') as fwd_out, open(output_fastq_rev, 'w+') as rev_out:
        with mp.Pool(processes=threads) as pool:
            # Process results as they complete
            for res in pool.imap_unordered(collapse_wrapper_paired, zip(bc_list, fwd_list, rev_list, min_list, flag_list)):
                if res is not None:
                    fwd_out.write(res[0])
                    rev_out.write(res[1])
    return None

def write_to_file_single(bc_list, fwd_list, output_fastq_fwd, min_count, align_flag, threads):
    min_list  = [min_count] * len(bc_list)
    flag_list = [align_flag] * len(bc_list)
    with open(output_fastq_fwd, 'w+') as fwd_out:
        with mp.Pool(processes=threads) as pool:
            # Process results as they complete
            for res in pool.imap_unordered(collapse_wrapper_single, zip(bc_list, fwd_list, min_list, flag_list)):
                if res is not None:
                    fwd_out.write(res)
    return None

def process_minibatch_paired(fwd_reads, fwd_qscores, rev_reads, rev_qscores, min_qscore, bc_start, bc_len, bc_end):
    # generally faster with mulitprocessing to hand each thread many reads at a time
    minibatch_fwd_map = {}
    minibatch_rev_map = {}
    for i in range(len(fwd_reads)):
        bc_int, fwd_int, rev_int = extract_read_paired(fwd_reads[i], fwd_qscores[i], rev_reads[i], rev_qscores[i], min_qscore, bc_start, bc_len, bc_end)
        if bc_int:
            if minibatch_fwd_map.get(bc_int) == None:
                minibatch_fwd_map[bc_int] = [fwd_int]
                minibatch_rev_map[bc_int] = [rev_int]
            else:
                minibatch_fwd_map[bc_int].append(fwd_int)
                minibatch_rev_map[bc_int].append(rev_int)
    return minibatch_fwd_map, minibatch_rev_map

def process_minibatch_single(fwd_reads, fwd_qscores, min_qscore, bc_start, bc_len, bc_end):
    # generally faster with mulitprocessing to hand each thread many reads at a time
    minibatch_fwd_map = {}
    for i in range(len(fwd_reads)):
        bc_int, fwd_int = extract_read_single(fwd_reads[i], fwd_qscores[i], min_qscore, bc_start, bc_len, bc_end)
        if bc_int:
            if minibatch_fwd_map.get(bc_int) == None:
                minibatch_fwd_map[bc_int] = [fwd_int]
            else:
                minibatch_fwd_map[bc_int].append(fwd_int)
    return minibatch_fwd_map

##########################################
#   Workhorse data-processing functions  #
##########################################

def extract_read_paired(fwd_read, fwd_qscore, rev_read, rev_qscore, min_qscore, bc_start, bc_len, bc_end):
    # Add quality filters here
    # check qscores
    for q in fwd_qscore[:-1]:
        if ord(q) < min_qscore:
            return None, None, None
    for q in rev_qscore[:-1]:
        if ord(q) < min_qscore:
            return None, None, None
    # Make sure you're not trying to grab a barcode from a read that's too short
    if len(fwd_read) <= bc_len:
        return None, None, None
    
    bc_int  = seq_to_twobit(fwd_read[bc_start:bc_end])
    fwd_int = seq_to_twobit(fwd_read[:-1])
    rev_int = seq_to_twobit(rev_read[:-1])
    return bc_int, fwd_int, rev_int

def extract_read_single(fwd_read, fwd_qscore, min_qscore, bc_start, bc_len, bc_end):
    # Add quality filters here
    # check qscores
    for q in fwd_qscore[:-1]:
        if ord(q) < min_qscore:
            return None, None

    # Make sure you're not trying to grab a barcode from a read that's too short
    if len(fwd_read) <= bc_len:
        return None, None
    
    bc_int  = seq_to_twobit(fwd_read[bc_start:bc_end])
    fwd_int = seq_to_twobit(fwd_read[:-1])
    return bc_int, fwd_int

def collapse_barcode_paired(twobit_bc, fwd_seqs, rev_seqs, min_count, align_flag):

    # check count and skip if below minimum
    count = len(fwd_seqs)
    if count < min_count:
        return None

    bc_seq = twobit_to_seq(twobit_bc)

    #if only one read, no need to align or build consensus
    if count == 1:
        fwd_seq = twobit_to_seq(fwd_seqs[0])
        rev_seq = twobit_to_seq(rev_seqs[0])
        fwd_out = "@fwd_consensus_bc_read;bc=%s;count=1\n%s\n+\n%s\n" % (
            bc_seq,
            fwd_seq,
            ''.join(["4" for q in range(len(fwd_seq))])
            )
        rev_out = "@rev_consensus_bc_read;bc=%s;count=1\n%s\n+\n%s\n" % (
            bc_seq,
            rev_seq,
            ''.join(["4" for q in range(len(rev_seq))])
            )
        return fwd_out, rev_out
    
    #else, do full processing
    # convert out of integer space
    fwd_seqs = [twobit_to_list(seq) for seq in fwd_seqs]
    rev_seqs = [twobit_to_list(seq) for seq in rev_seqs]
    
    # align if requested
    if align_flag:
        fwd_ali_seqs = [needleman_wunsch(fwd_seqs[0], seq) for seq in fwd_seqs]
        fwd_ali_seqs = multi_align_lists(fwd_ali_seqs)
        rev_ali_seqs = [needleman_wunsch(rev_seqs[0], seq) for seq in rev_seqs]
        rev_ali_seqs = multi_align_lists(rev_ali_seqs)
    else:
        fwd_ali_seqs = fwd_seqs
        rev_ali_seqs = rev_seqs
    
    # make consensus sequences
    fwd_cons_seq, fwd_cons_q = build_consensus(fwd_ali_seqs)
    rev_cons_seq, rev_cons_q = build_consensus(rev_ali_seqs)

    # format to fastq strings
    fwd_out = "@fwd_consensus_bc_read;bc=%s;count=%s\n%s\n+\n%s\n" % (
        bc_seq,
        count,
        ''.join([rev_bases[b] for b in fwd_cons_seq]),
        ''.join([chr(q + 33) for q in fwd_cons_q])
        )
    rev_out = "@rev_consensus_bc_read;bc=%s;count=%s\n%s\n+\n%s\n" % (
        bc_seq,
        count,
        ''.join([rev_bases[b] for b in rev_cons_seq]),
        ''.join([chr(q + 33) for q in rev_cons_q])
        )
    return fwd_out, rev_out

def collapse_barcode_single(twobit_bc, fwd_seqs, min_count, align_flag):

    # check count and skip if below minimum
    count = len(fwd_seqs)
    if count < min_count:
        return None

    bc_seq = twobit_to_seq(twobit_bc)

    #if only one read, no need to align or build consensus
    if count == 1:
        fwd_seq = twobit_to_seq(fwd_seqs[0])
        fwd_out = "@fwd_consensus_bc_read;bc=%s;count=1\n%s\n+\n%s\n" % (
            bc_seq,
            fwd_seq,
            ''.join(["4" for q in range(len(fwd_seq))])
            )
        return fwd_out

    #else, do full processing
    # convert out of integer space
    fwd_seqs = [twobit_to_list(seq) for seq in fwd_seqs]

    # align if requested
    if align_flag:
        fwd_ali_seqs = [needleman_wunsch(fwd_seqs[0], seq) for seq in fwd_seqs]
        fwd_ali_seqs = multi_align_lists(fwd_ali_seqs)
    else:
        fwd_ali_seqs = fwd_seqs

    # make consensus sequences
    fwd_cons_seq, fwd_cons_q = build_consensus(fwd_ali_seqs)

    # format to fastq strings
    fwd_out = "@fwd_consensus_bc_read;bc=%s;count=%s\n%s\n+\n%s\n" % (
        bc_seq,
        count,
        ''.join([rev_bases[b] for b in fwd_cons_seq]),
        ''.join([chr(q + 33) for q in fwd_cons_q])
        )

    return fwd_out

#################################
#   Utility / Helper Functions  #
#################################

def open_fastq_file(filename):
    """Opens a FASTQ file, handling gzip compression if necessary."""
    if filename.endswith('.fastq.gz'):
        return gzip.open(filename, 'rt')
    elif filename.endswith('.fastq'):
        return open(filename, 'r')
    else:
        raise ValueError("Unsupported file format. File must end with .fastq or .fastq.gz")

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
    seq_list = seq_list[::-1]
    return seq_list[:-1] #remove capping "C"

def multi_align_lists(ali_lists, max_iter = 10000):

    # Initialize output lists and reverse the input lists
    ali_lists = [list(reversed(seq)) for seq in ali_lists]  # Reverse only once
    out_lists = [[] for _ in ali_lists]
    
    # iterate through lists
    for _ in range(max_iter):
        if all(not seq for seq in ali_lists):  # Stop when all sequences are empty
            break

        # check if any sequence contains inserted base
        insertion_flag = any(seq and seq[-1] > 4 for seq in ali_lists)
        for i, seq in enumerate(ali_lists):
            if seq and (seq[-1] > 4 if insertion_flag else True):
                out_lists[i].append(seq.pop())
            else:
                out_lists[i].append(0)  # Add gap for sequences without insertions

    return out_lists

def collapse_wrapper_paired(args):
    return collapse_barcode_paired(*args)

def collapse_wrapper_single(args):
    return collapse_barcode_single(*args)

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

if __name__ == "__main__":
    main()

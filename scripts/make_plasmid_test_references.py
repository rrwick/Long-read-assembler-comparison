#!/usr/bin/env python3
"""
This script generates reference sequences which contain plasmids of various
number, size and depth. These can then be used as Badread input genomes for
the purpose of testing assemblers on plasmid-rich genomes.

This script could work with any input genomes, but I had the Orlek 2017
collection in mind when I wrote it:
https://www.sciencedirect.com/science/article/pii/S2352340917301567
"""

import argparse
import gzip
import itertools
import pathlib
import random
import sys


def get_arguments():
    parser = argparse.ArgumentParser(description='Make plasmid-containing reference sequences for assembler tests')

    # Options
    parser.add_argument('--plasmid_count', type=int, default=1000,
                        help='Total number of plasmids')
    parser.add_argument('--max_plasmids_per_genome', type=int, default=9,
                        help='Reference genomes with up to this many reference plasmids will be made')

    # Positional (required) arguments
    parser.add_argument('chromosome', type=str,
                        help='FASTA sequence of chromosome')
    parser.add_argument('all_plasmids', type=str,
                        help='FASTA containing all plasmids (e.g. nucleotideseq.fa from Orlek 2017)')
    parser.add_argument('out_dir', type=str,
                        help='Output directory')
    return parser.parse_args()


def main():
    args = get_arguments()
    chromosome_seq = load_chromosome_seq(args.chromosome)
    plasmid_seqs = load_fasta(args.all_plasmids)
    plasmid_seqs = choose_plasmids(plasmid_seqs, args.plasmid_count)

    i = 0
    file_counter = itertools.count(start=1)
    while i < len(plasmid_seqs):
        for group_count in range(1, args.max_plasmids_per_genome + 1):
            reference_sequences = [('chromosome', chromosome_seq, 1.0)]
            print('chromosome')
            for j in range(group_count):
                name, seq, seq_len = plasmid_seqs[i]
                depth = 10 ** random.uniform(-2.0, 1.0)
                reference_sequences.append((name, seq, depth))
                print(f'plasmid {j+1}: {name} length={seq_len} depth={depth:.3f}x')
                i += 1
                if i == len(plasmid_seqs):
                    break
            write_reference_to_file(reference_sequences, args.out_dir, file_counter)
            print()
            if i == len(plasmid_seqs):
                break


def load_chromosome_seq(chromosome_fasta_filename):
    chromosome_seq = load_fasta(chromosome_fasta_filename)
    assert len(chromosome_seq) == 1
    return chromosome_seq[0][1]


def choose_plasmids(plasmid_seqs, plasmid_count):
    """
    This function takes in all plasmid sequences and chooses a subset with a relative even
    length distribution.
    """
    chosen_plasmids = set()
    while len(chosen_plasmids) < plasmid_count:
        for size_small in range(0, 250000, 500):
            size_big = size_small + 500
            plasmids_in_range = [(name, seq, seq_len) for name, seq, seq_len in plasmid_seqs
                                 if size_small <= seq_len < size_big and name not in chosen_plasmids]
            if plasmids_in_range:
                name, seq, seq_len = random.choice(plasmids_in_range)
                chosen_plasmids.add(name)
    chosen_plasmids = list(chosen_plasmids)
    random.shuffle(chosen_plasmids)
    chosen_plasmids = set(chosen_plasmids[:plasmid_count])

    plasmid_seqs = [(name, seq, seq_len) for name, seq, seq_len in plasmid_seqs
                  if name in chosen_plasmids]
    random.shuffle(plasmid_seqs)
    return plasmid_seqs


def write_reference_to_file(reference_sequences, out_dir, file_counter):
    pathlib.Path(out_dir).mkdir(exist_ok=True)
    filename = f'{next(file_counter):03d}.fasta'
    filepath = pathlib.Path(out_dir) / pathlib.Path(filename)
    print(filepath)
    with open(filepath, 'wt') as fasta:
        for name, seq, depth in reference_sequences:
            fasta.write(f'>{name} length={len(seq)} depth={depth:.3f}x circular=true\n')
            fasta.write(f'{seq}\n')


def load_fasta(filename):
    fasta_seqs = []
    with get_open_func(filename)(filename, 'rt') as fasta_file:
        name, sequence = '', []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    seq = ''.join(sequence)
                    fasta_seqs.append((name, seq, len(seq)))
                    sequence = []
                name = line[1:]
            else:
                sequence.append(line)
        if name:
            seq = ''.join(sequence)
            fasta_seqs.append((name, seq, len(seq)))
    return fasta_seqs


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('Error: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('Error: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


if __name__ == '__main__':
    main()

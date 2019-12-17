#!/usr/bin/env python3
"""
This script prepares a bunch of genomes for use as Badread references.

For each output genome it:
  * sets the chromosome to 1x depth
  * sets any plasmids to a random depth (the min/max of the range depending on the plasmid size)
  * sets all output replicons to circular

This script is licensed under a Creative Commons Attribution 4.0 International License. You should
have received a copy of the license along with this work. If not, see
<http://creativecommons.org/licenses/by/4.0/>.
"""

import glob
import gzip
import math
import os
import random
import sys


def main():
    input_assembly_dir = sys.argv[1]
    output_assembly_dir = sys.argv[2]

    random.seed(0)
    os.makedirs(output_assembly_dir, exist_ok=True)

    chromosome_seqs, plasmid_seqs = {}, {}
    real_assemblies = glob.glob(input_assembly_dir + '/*.fna')

    for a in real_assemblies:
        seqs = load_fasta(a)
        assembly_name = a.split('/')[-1][:15]
        print(assembly_name)
        filename = output_assembly_dir + '/' + assembly_name + '.fasta'
        with open(filename, 'wt') as out_file:
            for name, seq in seqs:
                length = len(seq)
                if length > 400000:  # chromosome
                    depth = 1.0
                else:
                    depth = get_plasmid_depth(length)
                out_file.write(f'>{name} length={length} depth={depth:.2f}x circular=true\n')
                out_file.write(seq)
                out_file.write('\n')


def get_plasmid_depth(plasmid_length):
    upper_limit = (200000 / plasmid_length) + 1.0
    lower_limit = (-5000 / (plasmid_length + 10000)) + 0.5
    upper_limit_log = math.log(upper_limit)
    lower_limit_log = math.log(lower_limit)
    return math.exp(random.uniform(lower_limit_log, upper_limit_log))


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


def load_fasta(filename):
    fasta_seqs = []
    with get_open_func(filename)(filename, 'rt') as fasta_file:
        name = ''
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    fasta_seqs.append((name.split()[0], ''.join(sequence)))
                    sequence = []
                name = line[1:]
                short_name = name.split()[0]
            else:
                sequence.append(line)
        if name:
            fasta_seqs.append((name.split()[0], ''.join(sequence)))
    return fasta_seqs


if __name__ == '__main__':
    main()

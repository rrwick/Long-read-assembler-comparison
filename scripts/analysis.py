#!/usr/bin/env python3
"""
This script takes the following arguments (in this order):
  * assembly filename
  * reference genome fasta
  * mode: 'genome' or 'replicon'
"""

import collections
import gzip
import pathlib
import re
import subprocess
import sys
import tempfile

CONTIGUITY_COMPLETENESS_THRESHOLD = 0.99


def main():
    assembly_filename = sys.argv[1]
    ref_filename = sys.argv[2]
    mode = sys.argv[3]

    if not pathlib.Path(assembly_filename).is_file():
        sys.exit()
    if not pathlib.Path(ref_filename).is_file():
        sys.exit()

    assert mode == 'genome' or mode == 'replicon'
    genome_name = get_genome_name(assembly_filename)
    assembler = get_assembler_name(assembly_filename)

    ref_seqs = load_fasta(ref_filename)
    contigs = load_fasta(assembly_filename)
    replicon_names = [x[0] for x in ref_seqs]
    contig_names = [x[0] for x in contigs]

    with tempfile.TemporaryDirectory() as temp_dir:
        tripled_ref_filename = triple_reference(ref_seqs, temp_dir)
        contig_to_ref_alignments = get_contig_to_ref_alignments(assembly_filename, tripled_ref_filename)

    completed_chromosome, completed_all = True, True
    replicon_lines = []
    for replicon_name in replicon_names:
        replicon_length = [len(x[1]) for x in ref_seqs if x[0] == replicon_name][0]
        contig_name, contiguity, identity = get_contiguity(replicon_name, replicon_length, contig_to_ref_alignments)
        if contiguity < CONTIGUITY_COMPLETENESS_THRESHOLD:
            completed_all = False
            if replicon_length > 400000:
                completed_chromosome = False
        replicon_lines.append('\t'.join([genome_name, f'{replicon_name}', f'{replicon_length}',
                                         assembler, f'{contig_name}', f'{contiguity}', f'{identity}']))

    if mode == 'replicon':
        replicon_lines = sorted(replicon_lines, key=lambda x: int(x.split('\t')[2]), reverse=True)
        for line in replicon_lines:
            print(line)

    if mode == 'genome':
        completed_chromosome_binary = 1 if completed_chromosome else 0
        completed_all_binary = 1 if completed_all else 0
        print(f'{genome_name}\t{assembler}\t{completed_chromosome_binary}\t{completed_all_binary}')


def get_genome_name(assembly_filename):
    assembly_name = assembly_filename.split('/')[-1].replace('.fasta.gz', '')
    return '_'.join(assembly_name.rsplit('_')[:-2])


def get_assembler_name(assembly_filename):
    assembly_name = assembly_filename.split('/')[-1].replace('.fasta.gz', '')
    return '_'.join(assembly_name.rsplit('_')[-2:])


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
    try:
        fasta_seqs = []
        with get_open_func(filename)(filename, 'rt') as fasta_file:
            name = ''
            sequence = ''
            for line in fasta_file:
                line = line.strip()
                if not line:
                    continue
                if line[0] == '>':  # Header line = start of new contig
                    if name:
                        contig_name = name.split()[0]
                        fasta_seqs.append((contig_name, sequence))
                        sequence = ''
                    name = line[1:]
                else:
                    sequence += line
            if name:
                contig_name = name.split()[0]
                fasta_seqs.append((contig_name, sequence))
        return fasta_seqs
    except FileNotFoundError:
        return []


def triple_reference(ref_seqs, temp_dir):
    tripled_ref_filename = temp_dir + '/tripled.fasta'
    with open(tripled_ref_filename, 'wt') as tripled_ref:
        for name, seq in ref_seqs:
            tripled_ref.write(f'>{name}\n')
            tripled_ref.write(f'{seq}{seq}{seq}\n')
    return tripled_ref_filename


def get_contig_to_ref_alignments(assembly_filename, tripled_ref_filename):
    p = subprocess.run(['minimap2', '-c', '-t', '8', '-r', '10000', '-g', '10000', '-x', 'asm20', tripled_ref_filename, assembly_filename], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    minimap2_output = p.stdout.decode()
    alignments = []
    for line in minimap2_output.splitlines():
        alignments.append(Alignment(line))
    return alignments


def get_contiguity(ref_seq_name, replicon_length, contig_to_ref_alignments):
    alignments = [a for a in contig_to_ref_alignments if a.ref_name == ref_seq_name]
    if not alignments:
        return '', 0.0, 0.0
    best_alignment = sorted(alignments, key=lambda a: (a.ref_align_length, 1.0 / (a.ref_start+1)))[-1]
    contig_name = best_alignment.query_name
    contiguity = best_alignment.ref_align_length / replicon_length

    # To provide an identity value, the alignment must cover either 100 kbp or 1/4 of the replicon.
    if best_alignment.ref_align_length > 100000 or best_alignment.ref_align_length > replicon_length / 4:
        identity = best_alignment.identity
    else:
        identity = 0.0

    return contig_name, contiguity, identity


class Alignment(object):

    def __init__(self, paf_line):
        line_parts = paf_line.strip().split('\t')
        if len(line_parts) < 11:
            sys.exit('Error: alignment file does not seem to be in PAF format')

        self.query_name = line_parts[0]
        self.query_length = int(line_parts[1])
        self.query_start = int(line_parts[2])
        self.query_end = int(line_parts[3])
        self.strand = line_parts[4]

        self.ref_name = line_parts[5]
        self.ref_length = int(line_parts[6])
        self.ref_start = int(line_parts[7])
        self.ref_end = int(line_parts[8])

        self.identity = int(line_parts[9]) / int(line_parts[10])

        self.ref_align_length = self.ref_end - self.ref_start
        self.query_align_length = self.query_end - self.query_start

    def __repr__(self):
        return f'{self.query_name}: {self.query_start}-{self.query_end}; {self.ref_name}: {self.ref_start}-{self.ref_end}'


if __name__ == '__main__':
    main()

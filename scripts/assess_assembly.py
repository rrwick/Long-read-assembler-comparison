#!/usr/bin/env python3
"""
This script compares a bacterial whole genome assembly to a reference sequence. It assumes that
all replicons are circular, the reference sequence is complete (i.e. one contig per replicon) and
that the reference sequences are cleanly circularised (i.e. no missing or extra bases at the
start/end of the sequence).

It can be run in three ways:
1) assess_assembly.py --mode genome assembly.fasta reference.fasta
   When run in 'genome' mode, the output will contain a single line of results that indicates
   whether the chromosome was completely assembled and whether all replicons were completely
   assembled.
2) assess_assembly.py --mode replicon assembly.fasta reference.fasta
   When run in 'replicon' mode, the output will contain a line of results for each replicon in the
   reference genome which indicates the contiguity, identity and maximum indel size for the
   replicon.
3) assess_assembly.py --mode alignment assembly.fasta reference.fasta
   When run in 'alignment' mode, the output will show the entire longest alignment between the
   assembly and each reference sequence.
"""

import argparse
import gzip
import pathlib
import re
import subprocess
import sys
import tempfile


def get_arguments():
    parser = argparse.ArgumentParser(description='Assess a bacterial genome assembly against a '
                                                 'reference genome sequence')

    parser.add_argument('assembly_filename', type=str,
                        help='FASTA file to be assessed (can be gzipped)')
    parser.add_argument('ref_filename', type=str,
                        help='reference FASTA file (can be gzipped)')
    parser.add_argument('--mode', type=str, required=True,
                        help='Script mode: genome, replicon or alignment')

    parser.add_argument('--header', action='store_true',
                        help='Include a header line in the output')
    parser.add_argument('--threads', type=int, default=4,
                        help='Number of minimap2 alignment threads')
    parser.add_argument('--plasmid_threshold', type=int, default=400000,
                        help='Replicons less than or equal to this length are assumed to be '
                             'plasmids')
    parser.add_argument('--completeness_threshold', type=float, default=0.99,
                        help='Replicons with this much or more contiguity are considered complete')

    args = parser.parse_args()
    return args


def main():
    args = get_arguments()

    if not pathlib.Path(args.assembly_filename).is_file():
        sys.exit()
    if not pathlib.Path(args.ref_filename).is_file():
        sys.exit()

    assert args.mode == 'genome' or args.mode == 'replicon' or args.mode == 'alignment'
    assembly_name = pathlib.Path(args.assembly_filename).name

    ref_seqs = load_fasta(args.ref_filename)
    contigs = load_fasta(args.assembly_filename)
    replicon_names = [x[0] for x in ref_seqs]

    contig_to_ref_alignments = get_contig_to_ref_alignments(args.assembly_filename,
                                                            ref_seqs, args.threads)

    completed_chromosome, completed_all = True, True
    replicon_lines = []
    for replicon_name in replicon_names:
        replicon_length = [len(x[1]) for x in ref_seqs if x[0] == replicon_name][0]
        contig_name, contiguity, identity, max_indel = \
            get_contiguity(replicon_name, replicon_length, contig_to_ref_alignments, ref_seqs,
                           contigs, args.mode)
        if contiguity < args.completeness_threshold:
            completed_all = False
            if replicon_length > args.plasmid_threshold:
                completed_chromosome = False
        replicon_lines.append('\t'.join([assembly_name, f'{replicon_name}', f'{replicon_length}',
                                         f'{contig_name}', f'{contiguity}', f'{identity}',
                                         f'{max_indel}']))

    if args.mode == 'replicon':
        if args.header:
            print_replicon_header()
        replicon_lines = sorted(replicon_lines, key=lambda x: int(x.split('\t')[2]), reverse=True)
        for line in replicon_lines:
            print(line)

    elif args.mode == 'genome':
        if args.header:
            print_genome_header()
        completed_chromosome_binary = 1 if completed_chromosome else 0
        completed_all_binary = 1 if completed_all else 0
        print(f'{assembly_name}\t{completed_chromosome_binary}\t{completed_all_binary}')


def print_replicon_header():
    print('\t'.join(['assembly', 'replicon_name', 'length', 'contig_name', 'contiguity',
                     'identity', 'max_indel']))


def print_genome_header():
    print('\t'.join(['assembly', 'complete_chromosome', 'complete_all']))


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


def get_contig_to_ref_alignments(assembly_filename, ref_seqs, threads):
    with tempfile.TemporaryDirectory() as temp_dir:
        tripled_ref_filename = triple_reference(ref_seqs, temp_dir)
        p = subprocess.run(['minimap2', '-c', '-t', str(threads), '-x', 'asm20',
                            '-r', '10000', '-g', '10000', '-z', '1000,500',
                            tripled_ref_filename, assembly_filename],
                           stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        minimap2_output = p.stdout.decode()
        alignments = []
        for line in minimap2_output.splitlines():
            alignments.append(Alignment(line))
        return alignments


def get_contiguity(ref_seq_name, replicon_length, contig_to_ref_alignments, ref_seqs, contigs, mode):
    alignments = [a for a in contig_to_ref_alignments if a.ref_name == ref_seq_name]
    if not alignments:
        return '', 0.0, 0.0, ''
    best_alignment = sorted(alignments, key=lambda a: (a.ref_align_length, 1.0 / (a.ref_start+1)))[-1]
    contig_name = best_alignment.query_name
    contiguity = best_alignment.ref_align_length / replicon_length

    # To provide an identity value, the alignment must cover either 100 kbp or 1/4 of the replicon.
    if best_alignment.ref_align_length > 100000 or best_alignment.ref_align_length > replicon_length / 4:
        identity = best_alignment.identity
        max_indel = best_alignment.get_max_indel()
    else:
        identity = 0.0
        max_indel = ''

    if mode == 'alignment':
        ref_seq = [s for s in ref_seqs if s[0] == best_alignment.ref_name][0][1]
        ref_seq = ref_seq + ref_seq + ref_seq
        contig_seq = [s for s in contigs if s[0] == best_alignment.query_name][0][1]
        best_alignment.print_alignment(ref_seq, contig_seq)

    return contig_name, contiguity, identity, max_indel


REV_COMP_DICT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
                 'D': 'H', 'H': 'D', 'N': 'N', 'r': 'y', 'y': 'r', 's': 's', 'w': 'w', 'k': 'm',
                 'm': 'k', 'b': 'v', 'v': 'b', 'd': 'h', 'h': 'd', 'n': 'n', '.': '.', '-': '-',
                 '?': '?'}


def complement_base(base):
    try:
        return REV_COMP_DICT[base]
    except KeyError:
        return 'N'


def reverse_complement(seq):
    return ''.join([complement_base(x) for x in seq][::-1])


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

        self.cigar = None
        for part in line_parts:
            if part.startswith('cg:Z:'):
                self.cigar = part[5:]
        if self.cigar is None:
            sys.exit('Error: no CIGAR string found')

    def __repr__(self):
        return f'{self.query_name}: {self.query_start}-{self.query_end}; {self.ref_name}: {self.ref_start}-{self.ref_end}'

    def get_max_indel(self):
        max_indel = 0
        cigar_parts = re.findall(r'\d+\w', self.cigar)
        for cigar_part in cigar_parts:
            num = int(cigar_part[:-1])
            letter = cigar_part[-1]
            if (letter == 'I' or letter == 'D') and num > max_indel:
                max_indel = num
        return max_indel

    def print_alignment(self, ref_seq, query_seq):
        ref_seq = ref_seq[self.ref_start:self.ref_end]
        query_seq = query_seq[self.query_start:self.query_end]
        if self.strand == '-':
            query_seq = reverse_complement(query_seq)

        cigar_parts = re.findall(r'\d+\w', self.cigar)
        ref_pos, query_pos = 0, 0
        aligned_ref_seq, aligned_query_seq = [], []
        for cigar_part in cigar_parts:
            num = int(cigar_part[:-1])
            letter = cigar_part[-1]
            if letter == 'M':
                aligned_ref_seq.append(ref_seq[ref_pos:ref_pos+num])
                aligned_query_seq.append(query_seq[query_pos:query_pos+num])
                ref_pos += num
                query_pos += num
            if letter == 'I':
                aligned_ref_seq.append('-' * num)
                aligned_query_seq.append(query_seq[query_pos:query_pos+num])
                query_pos += num
            if letter == 'D':
                aligned_ref_seq.append(ref_seq[ref_pos:ref_pos+num])
                aligned_query_seq.append('-' * num)
                ref_pos += num

        aligned_ref_seq = ''.join(aligned_ref_seq)
        aligned_query_seq = ''.join(aligned_query_seq)

        print(f'{self.query_name} vs {self.ref_name} ({self.strand}):')
        print(aligned_query_seq)
        print(aligned_ref_seq)


if __name__ == '__main__':
    main()

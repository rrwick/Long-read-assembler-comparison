#!/usr/bin/env python3
"""
This script takes the following arguments (in this order):
  * assembly filename
  * read filename
  * alignment filename (minimap2 paf of assembly to reference)
  * reference genome fasta

It outputs lots of information about the read set and the assembly. Run it
with no arguments to get the header line.
"""

import collections
import gzip
import re
import subprocess
import sys
import tempfile

columns = ['read_filename', 'assembler', 'name', 'length', 'depth', 'contig_match', 'contig_name']

def main():
    try:
        assembly_filename = sys.argv[1]
        read_filename = sys.argv[2]
        ref_filename = sys.argv[3]

    # If no arguments were given, just print the header line.
    except IndexError:
        print('\t'.join(columns))
        sys.exit(0)

    print(assembly_filename, file=sys.stderr)
    short_read_filename = '/'.join(read_filename.split('/')[-2:])
    assembler = get_assembler_name(assembly_filename)

    ref_seqs = load_fasta(ref_filename)
    ref_seq_names = [x[0] for x in ref_seqs]

    with tempfile.TemporaryDirectory() as temp_dir:
        tripled_ref_filename = triple_reference(ref_seqs, temp_dir)
        read_to_ref_alignments = get_read_to_ref_alignments(read_filename, tripled_ref_filename)
        contig_to_ref_alignments = get_contig_to_ref_alignments(assembly_filename, tripled_ref_filename)

    for ref_seq_name in ref_seq_names:
        print(f'    {ref_seq_name}', file=sys.stderr)
        plasmid_depth = get_plasmid_depth(ref_seq_name, read_to_ref_alignments)
        plasmid_length = [len(x[1]) for x in ref_seqs if x[0] == ref_seq_name][0]
        contig_match, contig_name = check_for_contig_match(ref_seq_name, contig_to_ref_alignments)
        result = [short_read_filename, assembler,
                  f'{ref_seq_name}', f'{plasmid_length}', f'{plasmid_depth:.7f}',
                  f'{contig_match}', f'{contig_name}']
        print('\t'.join(result))


def get_lowest_window_identity(cigar, window_size):
    lowest_window_id = float('inf')
    expanded_cigar = get_expanded_cigar(cigar)
    for i in range(0, len(expanded_cigar) - window_size):
        cigar_window = expanded_cigar[i:i+window_size]
        window_id = cigar_window.count('=') / window_size
        if window_id < lowest_window_id:
            lowest_window_id = window_id
    if lowest_window_id == float('inf'):
        return 0.0
    return lowest_window_id


def get_expanded_cigar(cigar):
    expanded_cigar = []
    cigar_parts = re.findall(r'\d+[IDX=]', cigar)
    for cigar_part in cigar_parts:
        num = int(cigar_part[:-1])
        letter = cigar_part[-1]
        expanded_cigar.append(letter * num)
    return ''.join(expanded_cigar)


def get_assembler_name(assembly_filename):
    assembly_name = assembly_filename.split('/')[-1].replace('.fasta.gz', '')
    return assembly_name.split('_', 1)[1]


def load_fasta(filename):
    fasta_seqs = []
    with open(filename, 'rt') as fasta_file:
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


def triple_reference(ref_seqs, temp_dir):
    tripled_ref_filename = temp_dir + '/tripled.fasta'
    with open(tripled_ref_filename, 'wt') as tripled_ref:
        for name, seq in ref_seqs:
            tripled_ref.write(f'>{name}\n')
            tripled_ref.write(f'{seq}{seq}{seq}\n')
    return tripled_ref_filename


def get_read_to_ref_alignments(read_filename, tripled_ref_filename):
    """
    Returns a list containing the best alignment for each read.
    """
    p = subprocess.run(['minimap2', '-t', '8', '-x', 'map-ont', tripled_ref_filename, read_filename], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    minimap2_output = p.stdout.decode()
    alignments = {}
    for line in minimap2_output.splitlines():
        a = Alignment(line)
        if a.query_name not in alignments or alignments[a.query_name].ref_align_length < a.ref_align_length:
            alignments[a.query_name] = a
    return list(alignments.values())


def get_contig_to_ref_alignments(assembly_filename, tripled_ref_filename):
    p = subprocess.run(['minimap2', '-c', '-t', '8', '-r', '10000', '-g', '10000', '-x', 'asm20', tripled_ref_filename, assembly_filename], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    minimap2_output = p.stdout.decode()
    alignments = []
    for line in minimap2_output.splitlines():
        alignments.append(Alignment(line))
    return alignments


def get_plasmid_depth(ref_seq_name, read_to_ref_alignments):
    plasmid_alignments = [a for a in read_to_ref_alignments if a.ref_name == ref_seq_name]
    if not plasmid_alignments:
        return 0.0
    plasmid_lengths = set(a.ref_length // 3 for a in plasmid_alignments)
    assert len(plasmid_lengths) == 1
    plasmid_length = list(plasmid_lengths)[0]
    return sum(a.ref_end - a.ref_start for a in plasmid_alignments) / plasmid_length


def check_for_contig_match(ref_seq_name, contig_to_ref_alignments):
    """
    Checks to see if there is a contig which meets the following conditions:
      1) has an alignment to the tripled plasmid which covers >95% of the plasmid
    """
    plasmid_alignments = [a for a in contig_to_ref_alignments if a.ref_name == ref_seq_name]

    # Find the contig with the most alignment to the plasmid.
    bases_per_contig = collections.defaultdict(int)
    for a in plasmid_alignments:
        bases_per_contig[a.query_name] += a.ref_align_length
    if not bases_per_contig:
        return False, ''
    contig = max(bases_per_contig, key=bases_per_contig.get)
    contig_alignments = [a for a in plasmid_alignments if a.query_name == contig]

    # Get the best single alignment.
    if not contig_alignments:
        return False, ''
    best_alignment = sorted(contig_alignments, key=lambda a: (a.ref_align_length, 1.0 / (a.ref_start+1)))[-1]

    # Make sure the alignment covers enough of the plasmid.
    plasmid_length = best_alignment.ref_length // 3
    plasmid_coverage = best_alignment.ref_align_length / plasmid_length
    if plasmid_coverage < 0.95:
        return False, ''

    # Make sure the alignments covers enough of the contig.
    contig_length = best_alignment.query_length
    covered_bases = set()
    for a in contig_alignments:
        for i in range(a.query_start, a.query_end):
            covered_bases.add(i)
    contig_coverage = len(covered_bases) / contig_length
    if plasmid_coverage < 0.95:
        return False, ''

    return True, best_alignment.query_name


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

        self.ref_align_length = self.ref_end - self.ref_start
        self.query_align_length = self.query_end - self.query_start

    def __repr__(self):
        return f'{self.query_name}: {self.query_start}-{self.query_end}; {self.ref_name}: {self.ref_start}-{self.ref_end}'


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
This script takes a minimap2 PAF as input and outputs the mean depth
for each reference contig.
"""

import collections
import sys


def main():
    paf_filename = sys.argv[1]

    with open(paf_filename, 'rt') as paf:
        alignments = [Alignment(line) for line in paf]

    reference_lengths = {}
    aligned_bases = collections.defaultdict(int)
    for a in alignments:
        if a.ref_name not in reference_lengths:
            reference_lengths[a.ref_name] = a.ref_length
        else:
            assert reference_lengths[a.ref_name] == a.ref_length
        aligned_bases[a.ref_name] += (a.ref_end - a.ref_start)

    reference_lengths = sorted(reference_lengths.items(), key=lambda x: x[1], reverse=True)
    for ref_name, ref_length in reference_lengths:
        depth = aligned_bases[ref_name] / ref_length
        print(f'{ref_name}\t{depth:.3f}')


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

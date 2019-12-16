#!/usr/bin/env python3
"""
This script replaces all ambiguous bases in the assembly with A.

It assumes the input FASTA file has a one-line-per-sequence format.
"""

import glob
import math
import os
import random
import sys


def main():
    assembly_filename = sys.argv[1]

    bases = set()
    with open(assembly_filename, 'rt') as assembly:
        for header in assembly:
            assert header.startswith('>')
            seq = next(assembly).strip()
            bases |= set(b for b in seq)
    canonical_bases = set(['A', 'C', 'G', 'T'])
    ambiguous_bases = bases - canonical_bases
    if not ambiguous_bases:
        print(f'No ambiguous bases in {assembly_filename}')
        sys.exit(0)

    print(f'Replacing ambiguous bases in {assembly_filename}: {", ".join(sorted(ambiguous_bases))}')
    lines = []
    with open(assembly_filename, 'rt') as assembly:
        for header in assembly:
            seq = next(assembly).strip()
            for b in ambiguous_bases:
                seq = seq.replace(b, 'A')
            lines.append(header)
            lines.append(seq + '\n')
    with open(assembly_filename, 'wt') as assembly:
        for line in lines:
            assembly.write(line)


if __name__ == '__main__':
    main()

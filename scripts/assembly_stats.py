#!/usr/bin/env python3
"""
This script takes the following arguments (in this order):
  * assembly filename
  * read filename
  * alignment filename (minimap2 paf of assembly to reference)
  * reference genome length

It outputs lots of information about the read set and the assembly. Run it
with no arguments to get the header line.
"""

import dateutil.parser
import gzip
import re
import sys

columns = ['read_filename', 'depth',
           'glitch_rate', 'glitch_size', 'glitch_skip',
           'chimera_rate', 'junk_rate', 'random_rate',
           'start_adapter_length', 'end_adapter_length', 
           'fragment_length_mean', 'fragment_length_stdev', 'fragment_length_n50',
           'read_identity_mean', 'read_identity_max', 'read_identity_stdev',
           'assembler',
           'contiguity', 'identity', 'lowest_window_identity', 'coverage',
           'contigs', 'size', 'n50',
           'minutes']

def main():
    try:
        assembly_filename = sys.argv[1]
        read_filename = sys.argv[2]
        read_log = read_filename.replace('.fastq.gz', '.log')
        paf_filename = sys.argv[3]
        ref_length = int(sys.argv[4])

    # If no arguments were given, just print the header line.
    except IndexError:
        print('\t'.join(columns))
        sys.exit(0)

    short_read_filename = '/'.join(read_filename.split('/')[-2:])
    depth = get_depth(read_log, ref_length)
    glitch_rate, glitch_size, glitch_skip = get_glitch_levels(read_log)
    chimera_rate, junk_rate, random_rate = get_other_problems(read_log)
    start_adapter_length, end_adapter_length = get_adapter_lengths(read_log)
    fragment_length_mean, fragment_length_stdev, fragment_length_n50 = get_fragment_length(read_log)
    read_identity_mean, read_identity_max, read_identity_stdev = get_read_identity(read_log)
    assembler = get_assembler_name(assembly_filename)
    contiguity, identity, lowest_window_identity, coverage = get_alignment_stats(paf_filename, ref_length)
    contigs, size, n50 = get_assembly_stats(assembly_filename, ref_length)
    minutes = get_assembly_time(assembly_filename)

    result = [short_read_filename, f'{depth:.3f}',
              f'{glitch_rate}', f'{glitch_size}', f'{glitch_skip}',
              f'{chimera_rate:.5f}', f'{junk_rate:.5f}', f'{random_rate:.5f}',
              f'{start_adapter_length}', f'{end_adapter_length}',
              f'{fragment_length_mean}', f'{fragment_length_stdev}', f'{fragment_length_n50}',
              f'{read_identity_mean}', f'{read_identity_max}', f'{read_identity_stdev}',
              assembler,
              f'{contiguity:.7f}', f'{identity:.7f}', f'{lowest_window_identity:.7f}', f'{coverage:.7f}', 
              f'{contigs}', f'{size:.7f}', f'{n50:.7f}',
              f'{minutes:.2f}']
    print('\t'.join(result))


def get_alignment_stats(paf_filename, ref_length):
    # Tracks the longest single alignment, in terms of the reference bases. This can
    # exceed 100% in the case of overlapping ends.
    longest_alignment = 0
    longest_alignment_id = 0
    longest_alignment_cigar = ''

    # It also tracks the fraction of reference bases hit at least once (ranges from 0-100%).
    covered_bases = set()

    with open(paf_filename) as paf:
        for line in paf:
            parts = line.strip().split('\t')
            start, end = int(parts[7]), int(parts[8])
            matching_bases, total_bases = int(parts[9]), int(parts[10])
            cigar = [x for x in parts if x.startswith('cg:Z:')][0][5:]
            if end - start > longest_alignment:
                longest_alignment = end - start
                longest_alignment_id = matching_bases / total_bases
                longest_alignment_cigar = cigar
            longest_alignment = max(longest_alignment, end - start)
            for i in range(start, end):
                covered_bases.add(i % ref_length)

    relative_longest_alignment = longest_alignment / ref_length
    coverage = len(covered_bases) / ref_length
    lowest_window_id = get_lowest_window_identity(longest_alignment_cigar, 1000)

    return relative_longest_alignment, longest_alignment_id, lowest_window_id, coverage


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


def get_depth(read_log, ref_length):
    last_progress_line = None
    with open(read_log, 'rt') as log:
        for line in log:
            if line.startswith('Simulating: '):
                last_progress_line = line
    last_progress_line = last_progress_line.replace(',', '')
    last_progress_line = last_progress_line.split('reads')[-1].split(' bp')[0].strip()
    total_read_size = int(last_progress_line)
    return total_read_size / ref_length


def get_glitch_levels(read_log):
    glitch_rate, glitch_size, glitch_skip = None, None, None
    with open(read_log, 'rt') as log:
        for line in log:
            if 'Reads will have no glitches' in line:
                return 0, 0, 0
            if 'rate (mean distance between glitches)' in line:
                glitch_rate = int(line.split('=')[-1].strip())
            if 'size (mean length of random sequence)' in line:
                glitch_size = float(line.split('=')[-1].strip())
            if 'skip (mean sequence lost per glitch)' in line:
                glitch_skip = float(line.split('=')[-1].strip())
    assert glitch_rate is not None
    assert glitch_size is not None
    assert glitch_skip is not None
    return glitch_rate, glitch_size, glitch_skip


def get_other_problems(read_log):
    chimera_rate, junk_rate, random_rate = 0.0, 0.0, 0.0
    with open(read_log, 'rt') as log:
        for line in log:
            if 'chimera join rate:' in line:
                chimera_rate = float(line.split(': ')[-1].split('%')[0].strip()) / 100
            if 'junk read rate:' in line:
                junk_rate = float(line.split(': ')[-1].split('%')[0].strip()) / 100
            if 'random read rate:' in line:
                random_rate = float(line.split(': ')[-1].split('%')[0].strip()) / 100
    return chimera_rate, junk_rate, random_rate


def get_adapter_lengths(read_log):
    start_adapter_length, end_adapter_length = None, None
    with open(read_log, 'rt') as log:
        for line in log:
            if 'Start adapter:' in line:
                if 'Start adapter: none' in line:
                    start_adapter_length = 0
                else:
                    line = next(log)
                    assert 'seq:' in line
                    line = line.replace(' (randomly generated)', '').strip()
                    line = line.split('seq: ')[-1]
                    start_adapter_length = len(line)
            if 'End adapter:' in line:
                if 'End adapter: none' in line:
                    end_adapter_length = 0
                else:
                    line = next(log)
                    assert 'seq:' in line
                    line = line.replace(' (randomly generated)', '').strip()
                    line = line.split('seq: ')[-1]
                    end_adapter_length = len(line)
    assert start_adapter_length is not None
    assert end_adapter_length is not None
    return start_adapter_length, end_adapter_length


def get_fragment_length(read_log):
    fragment_length_mean, fragment_length_stdev, fragment_length_n50 = None, None, None
    with open(read_log, 'rt') as log:
        for line in log:
            if 'Generating fragment lengths from a gamma distribution:' in line:
                mean_line = next(log)
                stdev_line = next(log)
                n50_line = next(log)
                fragment_length_mean = int(mean_line.split('=')[1].split('bp')[0].strip())
                fragment_length_stdev = int(stdev_line.split('=')[1].split('bp')[0].strip())
                fragment_length_n50 = int(n50_line.split('=')[1].split('bp')[0].strip())
    assert fragment_length_mean is not None
    assert fragment_length_stdev is not None
    assert fragment_length_n50 is not None
    return fragment_length_mean, fragment_length_stdev, fragment_length_n50


def get_read_identity(read_log):
    read_identity_mean, read_identity_max, read_identity_stdev = None, None, None
    with open(read_log, 'rt') as log:
        for line in log:
            if 'Generating read identities from a beta distribution:' in line:
                mean_line = next(log)
                max_line = next(log)
                stdev_line = next(log)
                read_identity_mean = float(mean_line.split('=')[1].split('%')[0].strip()) / 100.0
                read_identity_max = float(max_line.split('=')[1].split('%')[0].strip()) / 100.0
                read_identity_stdev = float(stdev_line.split('=')[1].split('%')[0].strip()) / 100.0
    assert read_identity_mean is not None
    assert read_identity_max is not None
    assert read_identity_stdev is not None
    return read_identity_mean, read_identity_max, read_identity_stdev


def get_assembly_stats(assembly_filename, ref_length):
    contig_lengths = sorted(get_contig_lengths(assembly_filename), reverse=True)
    total_length = sum(contig_lengths)
    target_length = total_length * 0.5
    length_so_far = 0
    n50 = 0
    for contig_length in contig_lengths:
        length_so_far += contig_length
        if length_so_far >= target_length:
            n50 = contig_length
            break
    return len(contig_lengths), total_length / ref_length, n50 / ref_length


def get_assembly_time(assembly_filename):
    time_filename = assembly_filename.replace('.fasta.gz', '.time')
    times = []
    with open(time_filename, 'rt') as time_file:
        for line in time_file:
            times.append(dateutil.parser.parse(line.strip()))
    assert len(times) == 2
    elapsed_time = times[1] - times[0]
    return elapsed_time.seconds / 60.0


def get_contig_lengths(filename):
    lengths = []
    with gzip.open(filename, 'rt') as fasta_file:
        name = ''
        sequence = ''
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    lengths.append(len(sequence))
                    sequence = ''
                name = line[1:].split()[0]
            else:
                sequence += line
        if name:
            lengths.append(len(sequence))
    return lengths


if __name__ == '__main__':
    main()

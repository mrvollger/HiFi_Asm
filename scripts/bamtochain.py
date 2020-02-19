#!/usr/bin/env python3

# code originally from Pete Audano modified by Mitchell Vollger

"""
Parse CIGAR strings and generate a liftover chain file.
"""

import argparse
import gzip
import pandas as pd
import pysam
import sys

# CIGAR ID to string
AS_CIGAR_ID = {
    0: 'BAM_CMATCH',
    1: 'BAM_CINS',
    2: 'BAM_CDEL',
    3: 'BAM_CREF_SKIP',
    4: 'BAM_CSOFT_CLIP',
    5: 'BAM_CHARD_CLIP',
    6: 'BAM_CPAD',
    7: 'BAM_CEQUAL',
    8: 'BAM_CDIFF',
    9: 'BAM_CBACK',
    10: 'NM'
}

MATCH_SET = ['BAM_CMATCH', 'BAM_CEQUAL', 'BAM_CDIFF']

def encoder_plain(line):
    """
    String for plain files.

    :param line: Line to encode.

    :return: Encoded line (unaltered, plain files require no encoding).
    """
    return line


def encoder_bin(line):
    """
    String for binary files.

    :param line: Line to encode.

    :return: Encoded line.
    """

    return line.encode()


# Main
if __name__ == '__main__':

    # Parse arguments
    parser = argparse.ArgumentParser(description='Lift annotations to aligned contigs by CIGAR string')

    parser.add_argument('--in', '-i', dest='in_file', required=True,
                        help='Input alignment file')

    parser.add_argument('--out', '-o', dest='out_file', required=True,
                        help='Output chain file.')

    parser.add_argument('--fai', dest='fai_file', required=True,
                        help='Reference FAI file.')

    parser.add_argument('--read', dest='read_list', nargs='*',
                        help='List of aligned reads/contigs to include (all reads by default)')

    cmd_args = parser.parse_args()

    # Process arguments
    cmd_args.out_file = cmd_args.out_file.strip()
    cmd_args.in_file = cmd_args.in_file.strip()

    if cmd_args.read_list is not None:
        read_set = set(cmd_args.read_list)
    else:
        read_set = set()

    # Read chromosome sizes
    size_table = pd.read_csv(
        cmd_args.fai_file, header=None,
        names=('CHROM', 'SIZE', 'START', 'LINE_CHAR', 'LINE_BYTES'),
        usecols=('CHROM', 'SIZE'),
        index_col='CHROM',
        squeeze=True,
		sep="\t"
    )

    # Read input file
    record_count = 0

    # Get output function
    if cmd_args.out_file.lower().endswith('.gz'):
        out_func = gzip.open
        encoder = encoder_bin
    else:
        out_func = open
        encoder = encoder_plain

    # Parse BAM and write chain
    with out_func(cmd_args.out_file, 'w') as out_file:
        with pysam.AlignmentFile(cmd_args.in_file, 'r') as in_file:
            for record in in_file:

                record_count += 1
                sys.stderr.write("\rRecords: {}".format(record_count))

                # Process if read is in the read set or if the read set is empty
                if read_set and record.query_name not in read_set:
                    continue

                # Get length
                query_length = len(record.query_sequence)

                # Get CIGAR
                cigar_list = tuple((AS_CIGAR_ID[cigar_code], cigar_len) for cigar_code, cigar_len in record.cigar)

                if len(cigar_list) == 0:
                    print('Skipping record {}: No CIGAR'.format(record_count))

                # Count number of upstream clipped bases
                if cigar_list[0][0] == 'BAM_CSOFT_CLIP':
                    clip_up = cigar_list[0][1]

                    if cigar_list[1][0] not in MATCH_SET: #!= 'BAM_CMATCH'
                        raise RuntimeError(
                            'First CIGAR op after BAM_CSOFT_CLIP must be BAM_CMATCH: record={}, {}'.format(
                                record_count, record.cigarstring
                            )
                        )
                else:
                    if cigar_list[0][0] not in MATCH_SET: #!= 'BAM_CMATCH':
                        raise RuntimeError(
                            'First CIGAR op must be BAM_CSOFT_CLIP or BAM_CMATCH: record={}, {}'.format(
                                record_count, record.cigarstring
                            )
                        )

                    clip_up = 0

                # Count number of downstream clipped bases
                if cigar_list[-1][0] == 'BAM_CSOFT_CLIP':
                    clip_dn = cigar_list[-1][1]

                    if cigar_list[-2][0] not in MATCH_SET: #!= 'BAM_CMATCH':
                        raise RuntimeError(
                            'Last CIGAR op before BAM_CSOFT_CLIP must be BAM_CMATCH: record={}, {}'.format(
                                record_count, record.cigarstring
                            )
                        )
                else:
                    if cigar_list[-1][0] not in MATCH_SET: #!= 'BAM_CMATCH':
                        raise RuntimeError(
                            'Last CIGAR op must be BAM_CSOFT_CLIP or BAM_CMATCH: record={}, {}'.format(
                                record_count, record.cigarstring
                            )
                        )

                    clip_dn = 0

                # Write header
                out_file.write(encoder(
                    'chain 1000 {} {} + {} {} {} {} + {} {} {}\n'.format(
                        record.reference_name,
                        size_table[record.reference_name],
                        record.reference_start,
                        record.reference_start + record.reference_length,
                        record.query_name,
                        query_length,
                        clip_up,
                        query_length - clip_dn,
                        record_count
                    )
                ))

                # Write records
                block_record = None

                for cigar_op, cigar_len in cigar_list:

                    # Match
                    if cigar_op in MATCH_SET: #== 'BAM_CMATCH':

                        # Write previous block
                        if block_record is not None:
                            out_file.write(encoder('{}\t{}\t{}\n'.format(*block_record)))

                        block_record = [cigar_len, 0, 0]

                    elif cigar_op == 'BAM_CDEL':
                        block_record[1] += cigar_len

                    elif cigar_op == 'BAM_CINS':
                        block_record[2] += cigar_len

                    elif cigar_op in {'BAM_CSOFT_CLIP', 'BAM_CHARD_CLIP'}:
                        pass

                    else:
                        raise RuntimeError('Un-implemented CIGAR operation {}: record={}, {}'.format(
                            cigar_op, record_count, record.cigarstring
                        ))

                # Write last record
                if block_record is not None:
                    out_file.write(encoder('{}\n'.format(*block_record)))

                out_file.write(encoder('\n'))

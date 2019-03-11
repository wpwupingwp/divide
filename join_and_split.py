#!/usr/bin/python3

import argparse
import logging
from glob import glob
from sys import exit


# define logger
FMT = '%(asctime)s %(levelname)-8s %(message)s'
DATEFMT = '%I:%M:%S'
TEMP_LOG = 'Temp.log'
logging.basicConfig(format=FMT, datefmt=DATEFMT, level=logging.INFO,
                    handlers=[logging.StreamHandler(),
                              logging.FileHandler(TEMP_LOG)])
try:
    import coloredlogs
    coloredlogs.install(level=logging.INFO, fmt=FMT, datefmt=DATEFMT)
except ImportError:
    pass
log = logging.getLogger(__name__)


def join_fastq(arg):
    table = str.maketrans('ATCG', 'TAGC')
    split_text = arg.text
    split_quality = arg.quality

    def reverse_complement(seq):
        return seq.translate(table)[::-1]

    log.info(f'Forward file: {arg.left}')
    log.info(f'Reverse file: {arg.right}')
    out = '.'.join(arg.left.split('.')[:-1]) + '.merge'
    log.info(f'Merged file: {out}')
    log.warning('The program use "\\n" as line break.')
    left = open(arg.left, 'r')
    right = open(arg.right, 'r')
    output = open(out, 'w')
    # do not consider ambiguous base
    # assume fastq format is correct
    n = 0
    iter_zip = zip(left, right)
    for l, r in iter_zip:
        # use left id
        seq_id = l
        l, r = next(iter_zip)
        l_seq = l[:-1]
        r_seq = r[:-1]
        seq = f'{l_seq}{split_text}{reverse_complement(r_seq)}\n'
        third_line, _ = next(iter_zip)
        l, r = next(iter_zip)
        l_quality = l[:-1]
        r_quality = r[:-1]
        quality = f'{l_quality}{split_quality}{r_quality[::-1]}\n'
        output.write(seq_id)
        output.write(seq)
        output.write(third_line[:-1]+'\n')
        output.write(quality)
        n += 1
    left.close()
    right.close()
    output.close()
    log.info(f'{n} sequences were merged.')


def split_fastq(arg):
    table = str.maketrans('ATCG', 'TAGC')

    def reverse_complement(seq):
        return seq.translate(table)[::-1]

    files = list(glob(arg.merge))
    split_text = arg.text
    log.info(f'{len(files)} fastq files found.')
    for fastq in files:
        log.info(f'Split {fastq}')
        merged = open(fastq, 'r')
        _ = next(merged)
        l_seq, r_seq = next(merged)[:-1].split(split_text)
        if len(l_seq) == 0 or len(r_seq) == 0:
            log.critical('Empty sequence found. Please check the input or '
                         'split text!')
            merged.close()
            exit(-1)
        l_slice = slice(0, len(l_seq))
        r_slice = slice(len(l_seq)+len(split_text), None)
        log.info(f'Left slice: {repr(l_slice)}')
        log.info(f'Right slice: {repr(r_slice)}')
        merged.seek(0)
        left_f = fastq + '.f'
        right_f = fastq + '.r'
        left = open(left_f, 'w')
        right = open(right_f, 'w')
        log.info(f'Forward file: {left_f}')
        log.info(f'Reverse file: {right_f}')
        n = 0
        for line in merged:
            seq_id = line[:-1]
            raw_seq = next(merged)[:-1]
            l_seq = raw_seq[l_slice]
            r_seq = reverse_complement(raw_seq[r_slice])
            third_line = next(merged)[:-1]
            raw_quality = next(merged)[:-1]
            l_quality = raw_quality[l_slice]
            r_quality = raw_quality[r_slice][::-1]
            left.write(seq_id+'\n')
            left.write(l_seq+'\n')
            left.write(third_line+'\n')
            left.write(l_quality+'\n')
            right.write(seq_id+'\n')
            right.write(r_seq+'\n')
            right.write(third_line+'\n')
            right.write(r_quality+'\n')
            n += 1
        log.info(f'{n} records were splited.\n')


def parse_args():
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=main.__doc__)
    arg.add_argument('action', choices=('join', 'split'), help='action')
    arg.add_argument('-f', '--left', help='forward fastq file')
    arg.add_argument('-r', '--right', help='reverse fastq file')
    arg.add_argument('-m', '--merge', help='merged file/files')
    arg.add_argument('-t', '--text', default='JOINTEXT', help='join string')
    parsed = arg.parse_args()
    if parsed.action == 'join':
        if not all([parsed.left, parsed.right]):
            log.critical('Action "join" need two files as input!')
            arg.print_help()
            exit(-1)
    else:
        if parsed.merge is None:
            log.critical('Action "merge" need merge file/files as input!')
            arg.print_help()
            exit(-1)
    parsed.quality = 'A' * len(parsed.text)
    return parsed


def main():
    """
    Docstring.
    """
    arg = parse_args()
    log.info('Start')
    log.info(f'Action: {arg.action}')
    if arg.action == 'join':
        join_fastq(arg)
    else:
        split_fastq(arg)
    log.info('Done.')


if __name__ == '__main__':
    main()

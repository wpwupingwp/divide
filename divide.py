﻿#!/usr/bin/python3

import argparse
import os
import regex as re
from subprocess import run
from timeit import default_timer as timer
from Bio import SeqIO


def divide_by_barcode(merged, barcode_dict, arg):
    barcode_folder = os.path.join(arg.output, 'BARCODE')

    # edit it according to primer length
    barcode_info = {'total': 0,
                    'barcode_total': 0,
                    'head_barcode_mismatch': 0,
                    'head_barcode_mode_wrong': 0,
                    'head_tail_barcode_unequal': 0,
                    'tail_barcode_mismatch': 0,
                    'tail_barcode_mode_wrong': 0}
    #  parse_mode(mode):
    barcode_len, repeat = [int(i) for i in arg.mode.split('*')]
    barcode_full_len = barcode_len * repeat
    # analyze input files
    divided_files = set()
    handle_wrong = open(os.path.join(arg.output, 'barcode_wrong.fastq'), 'a')
    for record in SeqIO.parse(merged, 'fastq'):
        barcode_info['total'] += 1
        # ignore wrong barcode
        barcode_f = str(record.seq[:barcode_full_len])
        barcode_r = record.seq[-barcode_full_len:]
        # reverse complement for sequence end
        barcode_r = barcode_r.reverse_complement()
        barcode_r = str(barcode_r)
        barcode_split_f = list()
        barcode_split_r = list()
        if repeat == 1:
            barcode_split_f = [barcode_f, ]
            barcode_split_r = [barcode_r, ]
        else:
            for index in range(repeat):
                start = 0 + barcode_len * index
                barcode_split_f.append(barcode_f[start:(start+barcode_len)])
                barcode_split_r.append(barcode_r[start:(start+barcode_len)])
        # for default, only judge if barcode in 5' is right
        if barcode_f not in barcode_dict:
            barcode_info['head_barcode_mismatch'] += 1
            SeqIO.write(record, handle_wrong, 'fastq')
            continue
        # here use list.count
        if barcode_split_f.count(barcode_split_f[0]) != len(barcode_split_f):
            barcode_info['head_barcode_mode_wrong'] += 1
            SeqIO.write(record, handle_wrong, 'fastq')
            continue
        if arg.strict:
            if barcode_f != barcode_r:
                barcode_info['head_tail_barcode_unequal'] += 1
                SeqIO.write(record, handle_wrong, 'fastq')
                continue
            if barcode_r not in barcode_dict:
                barcode_info['tail_barcode_mismatch'] += 1
                SeqIO.write(record, handle_wrong, 'fastq')
                continue
            if barcode_split_r.count(barcode_split_r[0]) != len(
                    barcode_split_r):
                barcode_info['tail_barcode_mode_wrong'] += 1
                SeqIO.write(record, handle_wrong, 'fastq')
                continue
        name = barcode_info[barcode_f]
        output_file = os.path.join(barcode_folder, name+'.fastq')
        divided_files.add(output_file)
        with open(output_file, 'a') as handle:
            record.id = '{}-{}'.format(name, record.id)
            SeqIO.write(record, handle, 'fastq')
    handle_wrong.close()
    return barcode_info, divided_files, barcode_full_len


def divide_by_primer(divided_files, primer_info, arg, barcode_full_len):
    gene_folder = os.path.join(arg.output, 'GENE')
    barcode_gene_folder = os.path.join(arg.output, 'BARCODE-GENE')
    handle_wrong = open(os.path.join(arg.output, 'primer_wrong.fastq'), 'a')
    start = barcode_full_len + arg.adapter
    end = max([len(i) for i in primer_info.keys()]) + start
    divided_files = set()

    not_found = 0
    for fastq_file in divided_files:
        records = SeqIO.parse(fastq_file, 'fastq')
        for record in records:
            found = False
            string = str(record[start:end].seq)
            for pattern in primer_info:
                if re.search(pattern, string, re.BESTMATCH) is not None:
                    gene = primer_info[pattern]
                    found = True
                    break
            if not found:
                SeqIO.write(record, handle_wrong, 'fastq')
                not_found += 1
            barcode = record.id.split('-')[0]
            record.id = '{}-{}'.format(gene, record.id)
            handle_name = os.path.join(
                barcode_gene_folder, '{0}-{1}.fastq'.format(
                    barcode, gene))
            divided_files.add(handle_name)
            with open(handle_name, 'a') as handle:
                SeqIO.write(record, handle, 'fastq')
                handle_gene = open(os.path.join(gene_folder,
                                                '{}.fastq'.format(gene)), 'a')
                SeqIO.write(record, handle_gene, 'fastq')
    return divided_files, not_found


def flash(arg):
    # check flash
    files = arg.input
    output = arg.output
    if len(files) == 1:
        return files[0]
    elif len(files) == 2:
        for flash in ['flash2', 'flash']:
            check = run('{} --version'.format(flash), shell=True)
            if check == 0:
                run('{0} {1} {2} -d {3} -o out'.format(
                    flash, files[0], files[1], output), shell=True)
                return os.path.join(output, 'out.extendedFrags.fastq')
        raise Exception('FLASH not found and you give pair-end input!')
    else:
        raise Exception('Only support single or pair-end input!')


def get_barcode_info(arg):
    barcode_info = dict()
    with open(arg.barcode_file, 'r') as input_file:
        for line in input_file:
            if line.startswith(('barcode', 'Barcode')):
                continue
            line = line.strip().split(sep=',')
            barcode_info[line[0].upper()] = line[1]
    return barcode_info


def get_primer_info(arg):
    primer_dict = dict()
    with open(arg.primer_file, 'r') as input_file:
        for line in input_file:
            if line.startswith(('Gene', 'gene')):
                continue
            line = line.strip().split(sep=',')
            for i in line[1:]:
                pattern = re.compile(r'({}){{e<={}}}'.format(
                    i, arg.max_mismatch))
                primer_dict[pattern] = line[0]
    return primer_dict


def check_vsearch():
    # to be continued
    vsearch = 'vsearch'
    check = run('{} --version'.format(vsearch, shell=True))
    if check == 0:
        return vsearch
    else:
        return None


def process(paramter):
    (data, barcode, db_name, mode, strict, adapter, evalue, output) = paramter
    result = divide_run(data, barcode, db_name, mode, strict, adapter,
                        evalue, output)
    return result


def merge_dict(a, b):
    # merge b into a
    for key in b.keys():
        if key in a:
            a[key] += b[key]
        else:
            a[key] = b[key]
    return a


def parse_args():
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    arg.add_argument('-a', '--adapter', dest='adapter', default=14, type=int,
                     help='length of adapter')
    arg.add_argument('-b', dest='barcode_file',
                     help='csv file containing barcode info')
    arg.add_argument('-c', dest='threads', type=int, default=1,
                     help='CPU cores to use')
    arg.add_argument('-p', dest='primer_file',
                     help='csv file containing primer info')
    arg.add_argument('-j', '--join_by_n', action='store_true',
                     help=('if set, join sequences FLASH failed to merge by '
                           '"N"'*10))
    arg.add_argument('-s', '--strict', action='store_true',
                     help="if set, consider barcode on the 5' and 3'")
    arg.add_argument('-m', dest='mode', default='5*2',
                     help='''barcode mode, default value is 5*2, i.e.,
                        barcode with length 5 repeated 2 times''')
    arg.add_argument('input', nargs='+', help='input file, fastq format')
    arg.add_argument('-o', dest='output', default='Result', help='output path')
    return arg.parse_args()


def main():
    """
    Sequence likes this:
    [Barcode][Adapter][Primer][Sequence]
    Or:
    [Barcode][Adapter][Primer][Sequence][Barcode][Primer][Primer][Sequence]
    """
    start_time = timer()
    arg = parse_args()
    # reduce time cost by '.'
    # create folders
    barcode_folder = os.path.join(arg.output, 'BARCODE')
    gene_folder = os.path.join(arg.output, 'GENE')
    barcode_gene_folder = os.path.join(arg.output, 'BARCODE-GENE')
    try:
        os.mkdir(arg.output)
        os.mkdir(barcode_folder)
        os.mkdir(gene_folder)
        os.mkdir(barcode_gene_folder)
    except OSError:
        raise Exception('output exists, please use another name')

    barcode_dict = get_barcode_info(arg)
    primer_dict = get_primer_info(arg)

    # merge
    merged = flash(arg.input, arg.output)
    # divide barcode
    barcode_result, divided_files, barcode_full_len = divide_by_barcode(
        merged, barcode_dict, arg)
    result_files, primer_not_found = divide_by_primer(
        divided_files, primer_dict, arg, barcode_full_len)

    # write statistics
    barcode_info = os.path.join(arg.output, 'barcode_info.csv')
    with open(barcode_info, 'w') as handle:
        for record in barcode_result.items():
            handle.write('{0},{1} \n'.format(*record))
            handle.write('Primer not found, {}\n'.format(primer_not_found))

    end_time = timer()
    print('Finished with {0:.3f}s. You can find results in {1}.\n'.format(
        end_time-start_time, arg.output))


if __name__ == '__main__':
    main()

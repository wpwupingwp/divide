﻿#!/usr/bin/python3

import argparse
import os
from multiprocessing import Pool, cpu_count
from subprocess import call
from timeit import default_timer as timer
from core import divide_run


def flash(files, output):
    # check flash
    if len(files) == 1:
        return files[0]
    elif len(files) == 2:
        for flash in ['flash2', 'flash']:
            check = call('{} --version'.format(flash), shell=True)
            if check == 0:
                call('{0} {1} {2} -d {3} -o out'.format(
                    flash, files[0], files[1], output), shell=True)
                return os.path.join(output, 'out.extendedFrags.fastq')
        raise Exception('FLASH not found!')
    else:
        raise Exception('Only support single or pair-end input!')


def get_barcode_info(barcode_file):
    barcode_info = dict()
    with open(barcode_file, 'r') as input_file:
        for line in input_file:
            if line.startswith(('barcode', 'Barcode')):
                continue
            line = line.split(sep=',')
            barcode_info[line[0].upper()] = line[1].strip()
    return barcode_info


def get_primer_info(primer_file, output):
    primer = list()
    with open(primer_file, 'r') as input_file:
        for line in input_file:
            if line.startswith(('Gene', 'gene')):
                continue
            line = line.split(sep=',')
            primer.append([i.strip() for i in line])
    # join primer pairs
    primer_file = os.path.join(output, 'primer.fasta')
    with open(primer_file, 'w') as output:
        for line in primer:
            try:
                gene_name, forward, reverse = line
                output.write('>{0}\n{1}\n'.format(gene_name, reverse))
            except:
                gene_name, forward = line
            output.write('>{0}\n{1}\n'.format(gene_name, forward))
    return primer_file


def check_vsearch():
    vsearch = 'vsearch'
    check = call('{} --version'.format(vsearch, shell=True))
    if check == 0:
        return vsearch
    else:
        return None


def run(data):
    divide_run(data, barcode, db_name, arg.mode,
               arg.strict, arg.adapter, arg.evalue, arg.output)


def main():
    """
    Sequence likes this:
    [Barcode][Adapter][Primer][Sequence]
    Or:
    [Barcode][Adapter][Primer][Sequence][Barcode][Primer][Primer][Sequence]
    """
    start_time = timer()
    global arg
    arg = argparse.ArgumentParser()
    arg.add_argument('-a', '--adapter', dest='adapter', default=14, type=int,
                     help='length of adapter, typical 14 for AFLP')
    arg.add_argument('-b', dest='barcode_file',
                     help='csv file containing barcode info')
    arg.add_argument('-p', dest='primer_file',
                     help='csv file containing primer info')
    arg.add_argument('-e', dest='evalue', default=1e-5, type=float,
                     help='evalue for BLAST')
    arg.add_argument('-s', '--strict', action='store_true',
                     help="if set, consider barcode on the 5' and 3'")
    arg.add_argument('-m', dest='mode', default='5*2',
                     help='''barcode mode, default value is 5*2, i.e.,
                        barcode with length 5 repeated 2 times''')
    arg.add_argument('input', nargs='+', help='input file, fastq format')
    arg.add_argument('-o', dest='output', default='Result', help='output path')
    arg = arg.parse_args()
    # create folders
    barcode_folder = os.path.join(arg.output, 'BARCODE')
    gene_folder = os.path.join(arg.output, 'GENE')
    barcode_gene_folder = os.path.join(arg.output, 'BARCODE-GENE')
    try:
        os.mkdir(arg.output)
        os.mkdir(barcode_folder)
        os.mkdir(gene_folder)
        os.mkdir(barcode_gene_folder)
    except:
        raise Exception('output exists, please use another name')

    global barcode
    barcode = get_barcode_info(arg.barcode_file)
    global db_name
    primer_file = get_primer_info(arg.primer_file, arg.output)
    db_name = os.path.splitext(primer_file)[0]
    call('makeblastdb -in {0} -out {1} -dbtype nucl'.format(
        primer_file, db_name), shell=True)
    # split
    merged = flash(arg.input, arg.output)
    # parallel

    merged = [[i, merged] for i in range(cpu_count()-1)]
    pool = Pool(cpu_count()-1)
    result = pool.map(run, merged)
    pool.close()
    pool.join()

    print(result)
    end_time = timer()
    print('Finished with {0:.3f}s. You can find results in {1}.\n'.format(
        end_time-start_time, arg.output))


if __name__ == '__main__':
    main()

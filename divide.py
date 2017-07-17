#!/usr/bin/python3

import argparse
import os
from multiprocessing import Pool, cpu_count
from subprocess import call
from timeit import default_timer as timer
from core import divide_run
from glob import glob


def flash(files, output):
    # check flash
    if len(files) == 1:
        return False, files[0]
    elif len(files) == 2:
        for flash in ['flash2', 'flash']:
            check = call('{} --version'.format(flash), shell=True)
            if check == 0:
                call('{0} {1} {2} -d {3} -o out'.format(
                    flash, files[0], files[1], output), shell=True)
                return True, os.path.join(output, 'out.extendedFrags.fastq')
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


def run(paramter):
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
    arg = argparse.ArgumentParser()
    arg.add_argument('-a', '--adapter', dest='adapter', default=14, type=int,
                     help='length of adapter, typical 14 for AFLP')
    arg.add_argument('-b', dest='barcode_file',
                     help='csv file containing barcode info')
    arg.add_argument('-c', dest='core_number', type=int,
                     default=(cpu_count()-1), help='CPU cores to use')
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
    cores = arg.core_number
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

    barcode = get_barcode_info(arg.barcode_file)
    primer_file = get_primer_info(arg.primer_file, arg.output)
    db_name = os.path.splitext(primer_file)[0]
    call('makeblastdb -in {0} -out {1} -dbtype nucl'.format(
        primer_file, db_name), shell=True)
    # split
    is_merge, merged = flash(arg.input, arg.output)
    if is_merge:
        files = ['{}.{}'.format(merged, i) for i in range(cores)]
    else:
        files = ['{}.{}'.format(os.path.join(
            arg.output, merged), i) for i in range(cores)]
    with open(merged, 'r') as raw:
        handle = open(files[0], 'a')
        for n, line in enumerate(raw):
            # write 4k every time
            if n % 4000 == 0:
                handle = open(files[(n//4000) % cores], 'a')
            handle.write(line)
            handle.close()
    files = glob(os.path.join(arg.output, merged)+'.*')
    files = [(n, i) for n, i in enumerate(files)]
    # parallel
    merged = [(i, barcode, db_name, arg.mode, arg.strict, arg.adapter,
               arg.evalue, arg.output) for i in files]
    pool = Pool(cores)
    results = pool.map(run, merged)
    pool.close()
    pool.join()
    # merge result
    Barcode_info, Sample_info, Gene_info = results[0]
    for result in results[1:]:
        Barcode_info = merge_dict(Barcode_info, result[0])
        Sample_info = merge_dict(Sample_info, result[1])
        Gene_info = merge_dict(Gene_info, result[2])
    # write statistics
    barcode_info = os.path.join(arg.output, 'barcode_info.csv')
    with open(barcode_info, 'w') as handle:
        for record in Barcode_info.items():
            handle.write('{0},{1} \n'.format(*record))
    sample_info = os.path.join(arg.output, 'sample_info.csv')
    with open(sample_info, 'w') as handle:
        for record in Sample_info.items():
            handle.write('{0},{1} \n'.format(*record))
    gene_info = os.path.join(arg.output, 'gene_info.csv')
    with open(gene_info, 'w') as handle:
        for record in Gene_info .items():
            handle.write('{0},{1} \n'.format(*record))

    end_time = timer()
    print('Finished with {0:.3f}s. You can find results in {1}.\n'.format(
        end_time-start_time, arg.output))


if __name__ == '__main__':
    main()

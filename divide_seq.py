﻿#!/usr/bin/python3

import argparse
import datetime
import os
import regex as re
from subprocess import run, DEVNULL
from timeit import default_timer as timer
from Bio import SeqIO
from Bio.Data.IUPACData import ambiguous_dna_values


def divide_by_barcode(merged, barcode_dict, arg):
    barcode_folder = os.path.join(arg.output, 'BARCODE')
    barcode_info = {'Head barcode mismatch': 0,
                    'Head barcode mode wrong': 0,
                    'Head tail barcode conflict': 0,
                    'Tail barcode mismatch': 0,
                    'Tail barcode mode_wrong': 0,
                    'Total barcode': 0}
    #  parse_mode(mode):
    barcode_len, repeat = [int(i) for i in arg.mode.split('*')]
    barcode_full_len = barcode_len * repeat
    # analyze input files
    divided_files = set()
    handle_wrong = open(os.path.join(arg.output, 'barcode_wrong.fastq'), 'a')
    for record in SeqIO.parse(merged, 'fastq'):
        barcode_info['Total barcode'] += 1
        # ignore wrong barcode
        barcode_f = str(record.seq[:barcode_full_len])
        barcode_r = record.seq[-barcode_full_len:]
        # reverse complement for sequence end
        # or not????
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
        # here use list.count
        if barcode_split_f.count(barcode_split_f[0]) != len(barcode_split_f):
            barcode_info['Head barcode mode wrong'] += 1
            SeqIO.write(record, handle_wrong, 'fastq')
            continue
        # for default, only judge if barcode in 5' is right
        if barcode_f not in barcode_dict:
            barcode_info['Head barcode mismatch'] += 1
            SeqIO.write(record, handle_wrong, 'fastq')
            continue
        if arg.strict:
            if barcode_split_r.count(barcode_split_r[0]) != len(
                    barcode_split_r):
                barcode_info['Tail barcode mode wrong'] += 1
                SeqIO.write(record, handle_wrong, 'fastq')
                continue
            if barcode_r not in barcode_dict:
                barcode_info['Tail barcode mismatch'] += 1
                SeqIO.write(record, handle_wrong, 'fastq')
                continue
            if barcode_dict[barcode_f] != barcode_dict[barcode_r]:
                barcode_info['Head tail barcode conflict'] += 1
                SeqIO.write(record, handle_wrong, 'fastq')
                continue
        # use forward barcode to classify samples
        name = barcode_dict[barcode_f]
        output_file = os.path.join(barcode_folder, name+'.fastq')
        divided_files.add(output_file)
        with open(output_file, 'a') as handle:
            record.id = '{}-{}'.format(name, record.id)
            SeqIO.write(record, handle, 'fastq')
    handle_wrong.close()
    return barcode_info, divided_files, barcode_full_len


def divide_by_primer(divided_files, primer_info, arg, barcode_len, primer_len):
    gene_folder = os.path.join(arg.output, 'GENE')
    barcode_gene_folder = os.path.join(arg.output, 'BARCODE-GENE')
    handle_wrong = open(os.path.join(arg.output, 'primer_wrong.fastq'), 'a')
    start = barcode_len + arg.adapter
    end = start + primer_len
    result_files = set()
    primer_result = {i: 0 for i in primer_info.values()}

    not_found = 0
    for fastq_file in divided_files:
        records = SeqIO.parse(fastq_file, 'fastq')
        for record in records:
            found = False
            string = str(record[start:end].seq)
            for pattern in primer_info:
                if re.search(pattern, string) is not None:
                    gene = primer_info[pattern]
                    found = True
                    break
            if not found:
                SeqIO.write(record, handle_wrong, 'fastq')
                not_found += 1
                continue
            barcode = record.id.split('-')[0]
            primer_result[gene] += 1
            record.id = '-'.join([gene, record.id])
            if not arg.no_barcode:
                handle_name = os.path.join(
                    barcode_gene_folder, '{}-{}.fastq'.format(barcode, gene))
            else:
                handle_name = os.path.join(barcode_gene_folder, gene+'.fastq')
            result_files.add(handle_name)
            with open(handle_name, 'a') as handle:
                handle_gene = open(os.path.join(
                    gene_folder, '{}.fastq'.format(gene)), 'a')
                # cut barcode/adapter/primer
                if arg.cut:
                    SeqIO.write(record[end:-end],
                                handle_gene, 'fastq')
                    SeqIO.write(record[end:-end],
                                handle, 'fastq')
                else:
                    SeqIO.write(record, handle, 'fastq')
                    SeqIO.write(record, handle_gene, 'fastq')
    return primer_result, result_files, not_found


def flash(arg):
    # check flash
    files = arg.input
    output = arg.output
    if len(files) == 1:
        return files[0]
    elif len(files) == 2:
        check = run('flash --version', shell=True, stdout=DEVNULL)
        if check.returncode == 0:
            run('flash {} {} -d {} -o out'.format(
                files[0], files[1], output), shell=True)
            return os.path.join(output, 'out.extendedFrags.fastq')
        else:
            raise Exception('FLASH not found and you give pair-end input!')
    else:
        raise Exception('Only support single or pair-end input!')


def get_barcode_info(arg):
    barcode_info = dict()
    with open(arg.barcode_file, 'r') as input_file:
        for line in input_file:
            if line.startswith(('sample', 'Sample')):
                continue
            line = line.strip().split(sep=',')
            for i in line[1:]:
                barcode_info[i.upper()] = line[0]
    return barcode_info


def expand(seq):
    # nt_search in biopython do re.search
    # since here only need expanded seq, rewrite it
    result = ''
    for i in seq:
        try:
            value = ambiguous_dna_values[i]
        except KeyError:
            tprint('Illegal IUPAC ambiguous base "{}" found!.'.format(i))
            value = 'N'
        if len(value) == 1:
            result += value
        else:
            result += '[{}]'.format(value)
    return result


def get_primer_info(arg):
    primer_dict = dict()
    seq = list()
    with open(arg.primer_file, 'r') as input_file:
        for line in input_file:
            if line.startswith(('Gene', 'gene')):
                continue
            line = line.strip().split(sep=',')
            for i in line[1:]:
                seq.append(i)
                expand_i = expand(i.upper())
                pattern = re.compile(r'({}){{e<={}}}'.format(
                    expand_i, arg.max_mismatch), re.BESTMATCH)
                primer_dict[pattern] = line[0]
    # max primer len
    primer_len = max([len(i) for i in seq])
    return primer_dict, primer_len


def vsearch(fasta, arg):
    suffix1 = '.all_consensus'
    suffix2 = '.bigsize'
    command = ('vsearch --cluster_size {} --id {} --strand {} --sizeout '
               '--consout tmp.fasta --quiet').format(fasta, arg.id, arg.strand)
    run(command, shell=True)
    command2 = 'vsearch --sortbysize tmp.fasta --quiet --output ' '{}'.format(
        fasta+suffix1)
    run(command2, shell=True)
    command3 = ('vsearch --sortbysize {} --minsize {} --output {} '
                '--quiet'.format(fasta+suffix1, arg.minsize,
                                 fasta+suffix2))
    if arg.topn is not None:
        command3 += ' --topn {}'.format(arg.topn)
    run(command3, shell=True)


def parse_args():
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    arg.add_argument('input', nargs='+', help='input file, fastq format')
    arg.add_argument('-a', '--adapter', dest='adapter', default=14, type=int,
                     help='length of adapter')
    arg.add_argument('-b', dest='barcode_file',
                     help='csv file containing barcode info')
    arg.add_argument('-c', '--cut', action='store_true',
                     help="cut barcode in 5' and 3'")
    arg.add_argument('-d', dest='mode', default='5*2',
                     help='''barcode mode, default value is 5*2, i.e.,
                        barcode with length 5 repeated 2 times''')
    arg.add_argument('-m', dest='max_mismatch', type=int, default=3,
                     help='maximum mismatch in primer')
    arg.add_argument('-no_barcode', action='store_true',
                     help='skip dividing barcode')
    arg.add_argument('-p', dest='primer_file',
                     help='csv file containing primer info')
    arg.add_argument('-s', '--strict', action='store_true',
                     help="if set, consider barcode on the 5' and 3'")
    arg.add_argument('-o', dest='output', help='output path')

    vsearch = arg.add_argument_group('vsearch options')
    vsearch.add_argument('-no_vsearch', action='store_true',
                         help='skip vsearch')
    vsearch.add_argument('-consout', help='output file name')
    vsearch.add_argument('-id', type=float, default=0.97,
                         help='identity threshold')
    vsearch.add_argument('-minsize', type=int, default=5,
                         help='minimum abundance')
    vsearch.add_argument('-strand', choices=('plus', 'both'), default='both',
                         help='strand that cluster used,  plus or both')
    vsearch.add_argument('-topn', type=int, help='only keep best n seqs')
    # arg.print_help()
    return arg.parse_args()


def tprint(string):
    s = '{}\t{}'.format(datetime.datetime.now().time(), string)
    print(s, flush=True)


def main():
    """
    Sequence likes this:
    [Barcode][Adapter][Primer][Sequence]
    Or:
    [Barcode][Adapter][Primer][Sequence][Barcode][Primer][Adapter][Sequence]
    """
    start_time = timer()
    arg = parse_args()
    if arg.no_barcode:
        arg.adapter = 0
        arg.mode = None
    # create folders
    if arg.output is None:
        arg.output = arg.input[0].replace('.fastq', '')
    barcode_folder = os.path.join(arg.output, 'BARCODE')
    gene_folder = os.path.join(arg.output, 'GENE')
    barcode_gene_folder = os.path.join(arg.output, 'BARCODE-GENE')
    try:
        os.mkdir(arg.output)
        os.mkdir(barcode_folder)
        os.mkdir(gene_folder)
        os.mkdir(barcode_gene_folder)
    except OSError:
        print('Output exists, please use another name')
        raise
    # merge
    tprint('Divide start')
    tprint('Merging input ...')
    merged = flash(arg)
    tprint('Merge done')
    # divide barcode
    tprint('Dividing barcode ...')
    if not arg.no_barcode:
        barcode_dict = get_barcode_info(arg)
        barcode_result, divided_files, barcode_len = divide_by_barcode(
            merged, barcode_dict, arg)
        tprint('Dividing barcode finished.')
    else:
        barcode_result = {'Total barcode': 0}
        divided_files = [merged, ]
        barcode_len = 0
    # divide primer
    tprint('Dividing primer ...')
    primer_dict, primer_len = get_primer_info(arg)
    primer_result, result_files, primer_not_found = divide_by_primer(
        divided_files, primer_dict, arg, barcode_len, primer_len)
    tprint('Dividing genes finished.')
    # write statistics
    good_barcode = barcode_result['Total barcode'] * 2 - sum(
        barcode_result.values())
    barcode_info = os.path.join(arg.output, 'barcode_info.csv')
    with open(barcode_info, 'w') as handle:
        for record in barcode_result.items():
            handle.write('{0},{1}\n'.format(*record))
        handle.write('Good barcode,{}\n'.format(good_barcode))
    primer_info = os.path.join(arg.output, 'primer_info.csv')
    with open(primer_info, 'w') as handle:
        handle.write('Primer not found,{}\n'.format(primer_not_found))
        for record in primer_result.items():
            handle.write('{0},{1} \n'.format(*record))
    # vsearch
    tprint('Start vsearch.')
    if not arg.no_vsearch:
        check = run('vsearch --version', shell=True, stdout=DEVNULL,
                    stderr=DEVNULL)
        if check.returncode != 0:
            raise Exception('vsearch not found!')
        for i in result_files:
            print('>', end='', flush=True)
            vsearch(i, arg)
    print('')
    tprint('Done with vsearch.')
    end_time = timer()
    tprint('Divide done.')
    print('Finished with {0:.3f}s. You can find results in {1}.\n'.format(
        end_time-start_time, arg.output))


if __name__ == '__main__':
    main()

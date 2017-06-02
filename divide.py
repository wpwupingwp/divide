#!/usr/bin/python3

import argparse
import os
from Bio import SearchIO, SeqIO
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from collections import defaultdict
from multiprocessing import cpu_count
from subprocess import call
from timeit import default_timer as timer


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


def divide_barcode(merged, barcode, mode, strict,
                   adapter_length, output):
    # prepare input
    barcode_folder = os.path.join(output, 'BARCODE')
    os.mkdir(barcode_folder)

    # edit it according to primer length
    SEARCH_LEN = 25
    statistics = defaultdict(lambda: 0)
    #  parse_mode(mode):
    barcode_len, repeat = [int(i) for i in mode.split('*')]
    barcode_full_len = barcode_len * repeat
    skip = barcode_full_len + adapter_length
    # analyze input files
    divided_files = set()
    fastq_raw = SeqIO.parse(merged, 'fastq')
    handle_wrong = open(os.path.join(output, 'barcode_wrong.fastq'), 'w')
    head_file = os.path.join(output, 'head.fasta')
    handle_fasta = open(head_file, 'w')
    for record in fastq_raw:
        statistics['total'] += 1
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
            for index in range(repeat-1):
                start = 0 + barcode_len * index
                barcode_split_f.append(barcode_f[start:(start+barcode_len)])
                barcode_split_r.append(barcode_r[start:(start+barcode_len)])
        # for default, only judge if barcode in 5' is right
        if barcode_f not in barcode:
            statistics['head_barcode_mismatch'] += 1
            SeqIO.write(record, handle_wrong, 'fastq')
            continue
        # here use list.count
        if barcode_split_f.count(barcode_split_f[0]) != len(barcode_split_f):
            statistics['head_barcode_mode_wrong'] += 1
            SeqIO.write(record, handle_wrong, 'fastq')
            continue
        if strict:
            if barcode_f != barcode_r:
                statistics['head_tail_barcode_unequal'] += 1
                SeqIO.write(record, handle_wrong, 'fastq')
                continue
            if barcode_r not in barcode:
                statistics['tail_barcode_mismatch'] += 1
                SeqIO.write(record, handle_wrong, 'fastq')
                continue
            if barcode_split_r.count(barcode_split_r[0]) != len(
                    barcode_split_r):
                statistics['tail_barcode_mode_wrong'] += 1
                SeqIO.write(record, handle_wrong, 'fastq')
                continue
        name = barcode[barcode_f]
        output_file = os.path.join(barcode_folder, name+'.fastq')
        divided_files.add(output_file)
        with open(output_file, 'a') as handle:
            SeqIO.write(record, handle, 'fastq')
        handle_fasta.write('>{0}\n{1}\n'.format(
            record.description,
            record.seq[skip:(skip+SEARCH_LEN)]))
    handle_wrong.close()
    handle_fasta.close()
    return statistics, head_file, divided_files


def blast_and_parse(query_file, db_file, evalue, output):
    # build db
    db_name = os.path.splitext(db_file)[0]
    call('makeblastdb -in {0} -out {1} -dbtype nucl'.format(
        db_file, db_name), shell=True)
    blast_result_file = os.path.join(output, 'BlastResult.xml')
    cmd = nb(
        num_threads=cpu_count(),
        query=query_file,
        db=db_name,
        # Use blastn-short for primers.
        task='blastn-short',
        max_target_seqs=1,
        max_hsps=1,
        evalue=evalue,
        # xml output
        outfmt=5,
        out=blast_result_file
    )
    stdout, stderr = cmd()
    # parse
    blast_result = SearchIO.parse(blast_result_file, 'blast-xml')
    parse_result = dict()
    for record in blast_result:
        # skip empty blast result item
        if len(record) == 0:
            continue
        else:
            tophit = record[0]
        query_info = '{0} {1}'.format(
            tophit[0][0].query_id,
            tophit[0][0].query_description)
        hit_info = tophit[0][0].hit.id
        parse_result[query_info] = hit_info
    return parse_result


def divide_gene(blast_result, divided_files, output):
    gene_folder = os.path.join(output, 'GENE')
    barcode_gene_folder = os.path.join(output, 'BARCODE-GENE')
    os.mkdir(gene_folder)
    os.mkdir(barcode_gene_folder)
    # split and count
    sample_count = {i: 0 for i in divided_files}
    gene_count = defaultdict(lambda: 0)
    for fastq_file in divided_files:
        records = SeqIO.parse(fastq_file, 'fastq')
        for record in records:
            gene = record.description
            if gene in blast_result:
                sample_count[fastq_file] += 1
                gene_name = blast_result[gene]
                gene_count[gene_name] += 1
                barcode = os.path.splitext(fastq_file)[0]
                barcode = os.path.basename(barcode)
                record.id = '|'.join([gene_name, barcode, ''])
                handle_name = os.path.join(
                    barcode_gene_folder, '{0}-{1}.fastq'.format(
                        barcode, gene_name))
                with open(handle_name, 'a') as handle:
                    SeqIO.write(record, handle, 'fastq')
                    handle_gene = open(os.path.join(
                        gene_folder, '{}.fastq'.format(gene_name)), 'a')
                    SeqIO.write(record, handle_gene, 'fastq')
    return sample_count, gene_count


def main():
    """
    Sequence likes this:
    [Barcode][Adapter][Primer][Sequence]
    Or:
    [Barcode][Adapter][Primer][Sequence][Barcode][Primer][Primer][Sequence]
    """
    start_time = timer()
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
    try:
        os.mkdir(arg.output)
    except:
        raise Exception('output exists, please use another name')

    merged = flash(arg.input, arg.output)
    barcode = get_barcode_info(arg.barcode_file)
    primer_file = get_primer_info(arg.primer_file, arg.output)
    barcode_info, head_file, divided_files = divide_barcode(merged,
                                                            barcode,
                                                            arg.mode,
                                                            arg.strict,
                                                            arg.adapter,
                                                            arg.output)
    parse_result = blast_and_parse(head_file, primer_file,
                                   arg.evalue, arg.output)
    sample_info, gene_info = divide_gene(parse_result,
                                         divided_files, arg.output)
    with open(os.path.join(arg.output, 'barcode_info.csv'), 'w') as handle:
        for record in barcode_info.items():
            handle.write('{0},{1} \n'.format(*record))
    with open(os.path.join(arg.output, 'sample_info.csv'), 'w') as handle:
        for record in sample_info.items():
            handle.write('{0},{1} \n'.format(*record))
    with open(os.path.join(arg.output, 'gene_info.csv'), 'w') as handle:
        for record in gene_info .items():
            handle.write('{0},{1} \n'.format(*record))

    end_time = timer()
    print('Finished with {0:.3f}s. You can find results in {1}.\n'.format(
        end_time-start_time, arg.output))


if __name__ == '__main__':
    main()

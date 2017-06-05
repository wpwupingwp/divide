#!/usr/bin/python3

import os
from Bio import SearchIO, SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline as nb
from collections import defaultdict
from multiprocessing import cpu_count
from subprocess import call


def divide_run(data, barcode, primer_file, mode, strict,
               adapter_length, evalue, output):
    # prepare input
    thread_id, merged = data
    print(thread_id, merged)
    barcode_folder = os.path.join(output, 'BARCODE')

    # edit it according to primer length
    SEARCH_LEN = 25
    barcode_info = defaultdict(lambda: 0)
    #  parse_mode(mode):
    barcode_len, repeat = [int(i) for i in mode.split('*')]
    barcode_full_len = barcode_len * repeat
    skip = barcode_full_len + adapter_length
    # analyze input files
    divided_files = set()
    fastq_raw = SeqIO.parse(merged, 'fastq')
    handle_wrong = open(os.path.join(output, 'barcode_wrong.fastq'), 'a')
    head_file = os.path.join(output, 'head.fasta')
    handle_fasta = open(head_file, 'a')
    for record in fastq_raw:
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
            for index in range(repeat-1):
                start = 0 + barcode_len * index
                barcode_split_f.append(barcode_f[start:(start+barcode_len)])
                barcode_split_r.append(barcode_r[start:(start+barcode_len)])
        # for default, only judge if barcode in 5' is right
        if barcode_f not in barcode:
            barcode_info['head_barcode_mismatch'] += 1
            SeqIO.write(record, handle_wrong, 'fastq')
            continue
        # here use list.count
        if barcode_split_f.count(barcode_split_f[0]) != len(barcode_split_f):
            barcode_info['head_barcode_mode_wrong'] += 1
            SeqIO.write(record, handle_wrong, 'fastq')
            continue
        if strict:
            if barcode_f != barcode_r:
                barcode_info['head_tail_barcode_unequal'] += 1
                SeqIO.write(record, handle_wrong, 'fastq')
                continue
            if barcode_r not in barcode:
                barcode_info['tail_barcode_mismatch'] += 1
                SeqIO.write(record, handle_wrong, 'fastq')
                continue
            if barcode_split_r.count(barcode_split_r[0]) != len(
                    barcode_split_r):
                barcode_info['tail_barcode_mode_wrong'] += 1
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
    barcode_info_file = os.path.join(output, 'barcode_info.csv')
    with open(barcode_info_file, 'w') as stat_out:
        for record in barcode_info.items():
            stat_out.write('{0},{1} \n'.format(*record))
#     return statistics, head_file, divided_files


# def blast_and_parse(query_file, db_file, evalue, output):
    # build db
    db_name = os.path.splitext(primer_file)[0]
    call('makeblastdb -in {0} -out {1} -dbtype nucl'.format(
        primer_file, db_name), shell=True)
    blast_result_file = os.path.join(output,
                                     'BlastResult.xml.{}'.format(thread_id))
    cmd = nb(
        # num_threads=cpu_count(),
        query=head_file,
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


# def divide_gene(parse_result, divided_files, output):
    gene_folder = os.path.join(output, 'GENE')
    barcode_gene_folder = os.path.join(output, 'BARCODE-GENE')
    # split and count
    sample_count = {i: 0 for i in divided_files}
    gene_count = defaultdict(lambda: 0)
    for fastq_file in divided_files:
        records = SeqIO.parse(fastq_file, 'fastq')
        for record in records:
            gene = record.description
            if gene in parse_result:
                sample_count[fastq_file] += 1
                gene_name = parse_result[gene]
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
#    return sample_count, gene_count
    # write statistics
    sample_info = os.path.join(output, 'sample_info.csv')
    with open(sample_info, 'w') as handle:
        for record in sample_count.items():
            handle.write('{0},{1} \n'.format(*record))
    gene_info = os.path.join(output, 'gene_info.csv')
    with open(gene_info, 'w') as handle:
        for record in gene_count .items():
            handle.write('{0},{1} \n'.format(*record))

    return barcode_info, sample_info, gene_info

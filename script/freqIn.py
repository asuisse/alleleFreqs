from __future__ import division
from __future__ import print_function

import fnmatch
import json
import ntpath
import os
import sys
from collections import defaultdict
from optparse import OptionParser

import pandas as pd
import vcf


def find_normal(options):
    sample = ntpath.basename(options.variants).split("_")[0]
    config = pd.read_csv(options.config, delimiter="\t", index_col=False, na_filter=False, names=['sample', 'assay'])

    samples = config['sample'].tolist()
    it = iter(samples)

    s = {}
    for x in it:
        s[x] = next(it)

    if s[sample]:
        print("Tumour: %s" % sample)
        print("Normal: %s" % s[sample])
        return sample, s[sample]
    else:
        print("Cant find corresponding normal sample for %s" % sample)


def calculate_snp_frequency(options, tumour, normal, chrom, start, end, vcf_reader, supporting_records, opposing_records):
    total_freq, snp_count, support, oppose, total_af, ratio = 0, 0, 0, 0, 0, 0
    # print("Looking for snps in region: %s:%s-%s" % (chrom, start, end))
    mut_types = defaultdict(int)

    for record in vcf_reader.fetch(chrom, start, end):
        if record.genotype(normal)['DP'] > 20 and record.genotype(tumour)['GQ'] > 1:  # Don't filter on TUM DP to keep in DELS
            if (record.genotype(tumour)['GT'] == '0/0' and record.genotype(normal)['GT'] == '0/0') or (record.genotype(tumour)['GT'] == '1/1' and record.genotype(normal)['GT'] == '1/1'):
                continue
            if 'snp' not in record.INFO['TYPE'] or len(record.INFO['TYPE']) > 1:
                continue

            try:
                status = record.INFO['VT']
            except KeyError:
                status = 'germline'
                pass

            mut_types[status] += 1
            snp_count += 1
            accepted_genotypes = ['0/0', '0/1', '1/0', '1/1']

            if record.genotype(tumour)['GT'] not in accepted_genotypes or record.genotype(normal)['GT'] not in accepted_genotypes:
                continue

            # Throws error if any vals == 0
            taf = round((record.genotype(tumour)['AO'] / (record.genotype(tumour)['AO'] + record.genotype(tumour)['RO'])), 2)
            naf = round((record.genotype(normal)['AO'] / (record.genotype(normal)['AO'] + record.genotype(normal)['RO'])), 2)
            af_diff = abs(naf - taf)
            total_af += af_diff

            t_a_count = record.genotype(tumour)['AO']
            n_a_count = record.genotype(normal)['AO']

            if isinstance(t_a_count, list): t_a_count = t_a_count[0]  # Some records contain two values ... (?)
            if isinstance(n_a_count, list): n_a_count = n_a_count[0]  # Some records contain two values ... (?)

            t_freq = round((t_a_count / record.genotype(tumour)['DP']) * 100, 2)
            n_freq = round((n_a_count / record.genotype(normal)['DP']) * 100, 2)

            total_freq, support, supporting_records, oppose, opposing_records = is_shift(options, tumour, t_freq, normal, n_freq, total_freq, support, supporting_records, oppose, opposing_records, record)


    # if (snp_count / (end-start)) >= 0.005:
    av_freq = 0
    if snp_count > 0:
        av_freq = round(total_freq / snp_count, 2)
        print("Average across region: %s [ %s/%s ]" % (av_freq, total_freq, snp_count))
        print("SNVs in support of shift: %s [%s oppose]" % (support, oppose))

    sstring = "snp_su:%s" % support
    ostring = "snp_op:%s" % oppose

    if support > 2 or oppose > 2:
        ratio = round(((support + 0.01) / (oppose + 0.01)), 2)
        ratio_string = "snp_ratio:%s" % ratio
        nstring = '; '.join([sstring, ostring, ratio_string])
    else:
        nstring = '; '.join([sstring, ostring])

    if 'somatic' in mut_types:
        print("There's at least 1 somatic mutation in this region")
        ssnps = "somatic_snvs:%s" % mut_types['somatic']
        nstring = '; '.join([nstring, ssnps])

    # print(json.dumps(mut_types, indent=4, sort_keys=True))
    sites_surveyed = support + oppose
    return av_freq, nstring, supporting_records, sites_surveyed, ratio, opposing_records


def is_shift(options, tumour, t_freq, normal, n_freq, total_freq, support, supporting_records, oppose, opposing_records, record):
    difference = abs(t_freq - n_freq)

    if difference > 10:
        if options.event:
            print_event_details(difference, tumour, t_freq, normal, n_freq, record, True)
        total_freq += difference
        support += 1

        if options.write_vcf:
            # record = vcf.model._Record(CHROM=record.CHROM, POS=record.POS, ID='.', REF=record.REF, ALT=record.ALT, QUAL='.', FILTER='PASS', INFO={"TF": t_freq, "GT": }, FORMAT='.', sample_indexes=[], samples=None)
            supporting_records.append(record)  # vcf_writer.write_record(record)
    else:
        if options.event:
            print_event_details(difference, tumour, t_freq, normal, n_freq, record, False)
        oppose += 1
        if options.write_vcf:
            # record = vcf.model._Record(CHROM=record.CHROM, POS=record.POS, ID='.', REF=record.REF, ALT=record.ALT, QUAL='.', FILTER='PASS', INFO={"TF": t_freq, "GT": }, FORMAT='.', sample_indexes=[], samples=None)
            opposing_records.append(record)  # vcf_writer.write_record(record)

    return total_freq, support, supporting_records, oppose, opposing_records


def print_event_details(difference, tumour, t_freq, normal, n_freq, record, status):
    print("    ---- [ %s:%s ] ----" % (record.CHROM, record.POS))
    if status:
        print("    * --> This is a shift: %s " % difference)
        print("    * Tumour alt freq: %s" % t_freq)
        print("    * Normal alt freq: %s" % n_freq)
        print("      Tum depth: %s\tref: %s\talt: %s" % (record.genotype(tumour)['DP'], record.genotype(tumour)['RO'], record.genotype(tumour)['AO']))
        print("      Norm depth: %s\tref: %s\talt: %s" % (record.genotype(normal)['DP'], record.genotype(normal)['RO'], record.genotype(normal)['AO']))
    else:
        print("    --> This is NOT a shift: %s " % difference)
        print("      Tumour alt freq: %s" % t_freq)
        print("      Normal alt freq: %s" % n_freq)


def extract_vars(options):
    breakpoints = pd.read_csv(options.variants, delimiter="\t", index_col=False, na_filter=False)

    if options.event:
        breakpoints = breakpoints.loc[breakpoints['event'] == int(options.event)]

    tumour, normal = find_normal(options)

    # directory = '/Users/Nick_curie/Desktop/script_test/alleleFreqs/data'
    directory = '/Volumes/perso/Analysis/Analysis/Freebayes/vcf/'
    snps_file = '_'.join([tumour, 'snps_filt.vcf.gz'])
    snps = os.path.join(directory, snps_file)

    print("Reading snps from %s" % snps)

    vcf_reader = vcf.Reader(open(snps, 'r'))

    supporting_records = []
    opposing_records = []

    for index, row in breakpoints.iterrows():
        av_freq, sites_surveyed, ratio = 0, 0, 0
        nstring = None
        if row['type'] != "TRA" and row['chromosome1'] == row['chromosome2']:
            if (abs(row['log2(cnv)']) > 0.25 or (row['split_reads'] == '-' and row['disc_reads'] == '-')):
                print(("------------------------ Variant: %s %s [%s] ------------------------") % (row['event'], row['type'], row['position']))
                av_freq, nstring, supporting_records, sites_surveyed, ratio, opposing_records = calculate_snp_frequency(options, tumour, normal, row['chromosome1'], row['bp1'], row['bp2'], vcf_reader, supporting_records, opposing_records)

        length = row['bp2'] - row['bp1']

        breakpoints.loc[index, 'snp_freq'] = av_freq

        if nstring: write_notes(breakpoints, index, nstring)


        snp_support = do_snps_support(av_freq, length, sites_surveyed)

        if snp_support and sites_surveyed/length >= 0.0002 and row['status'] != '' and ratio >= 2:
            breakpoints.loc[index, 'status'] = 'T'

        elif not snp_support and sites_surveyed > 5 and sites_surveyed/length >= 0.0002 and row['status'] == '' and ('DUP' in row['type'] or 'DEL' in row['type']) and row['split_reads'] == '-' and row['disc_reads'] == '-':
            print("False positive. snp support: %s, sites surveyed: %s, length: %s" % (snp_support, sites_surveyed, length))
            breakpoints.loc[index, 'status'] = 'aF'

        print("Event: %s [%s], snp support: %s, sites surveyed: %s, length: %s, ratio: %s" % (row['event'], row['type'], snp_support, sites_surveyed, length, ratio))

    if not options.event:
        options.out_file = tumour + '_snv_added.txt'
        breakpoints.to_csv(options.out_file, sep="\t", index=False)

    if options.write_vcf:
        writeVCF(tumour, supporting_records, opposing_records, vcf_reader)


def write_notes(breakpoints, index, nstring):
    if breakpoints.loc[index, 'notes'] == '-' or breakpoints.loc[index, 'notes'] == '':
        breakpoints.loc[index, 'notes'] = nstring
    else:
        breakpoints.loc[index, 'notes'] = breakpoints.loc[index, 'notes'] + "; " + nstring


def do_snps_support(av_freq, length, sites_surveyed):
    if av_freq >= 10 and sites_surveyed >= 5:
        return True
    return False


def writeVCF(tumour, supporting_records, opposing_records, vcf_reader):

    sup_vcf_out = '_'.join([tumour, 'supporting_snps.vcf'])
    op_vcf_out = '_'.join([tumour, 'opposing_snps.vcf'])

    sup_vcf = vcf.Writer(open(sup_vcf_out, 'w'), vcf_reader)
    op_vcf = vcf.Writer(open(op_vcf_out, 'w'), vcf_reader)

    for record in supporting_records:
        sup_vcf.write_record(record)

    for record in opposing_records:
        op_vcf.write_record(record)


def get_args():
    parser = OptionParser()

    parser.add_option("-v", "--variants", dest="variants", action="store", help="svParser format file")
    parser.add_option("--event", dest="event", action="store", help="Inspect one event")
    parser.add_option("--config", dest="config", action="store", help="mapping for tumour/normal samples")
    parser.add_option("-o", "--out_file", dest="out_file", action="store", help="File to write annotated vars to")
    parser.add_option("--write_vcf", dest="write_vcf", action="store_true", help="VCF file to write supporting snvs to")
    parser.set_defaults(config='/Users/Nick_curie/Desktop/script_test/alleleFreqs/data/samples.tsv')

    options, args = parser.parse_args()

    if not options.variants:
        parser.print_help()
        print ()
        sys.exit("[!] Must provide a variants file. Exiting.")
    else:
        return options, args


def main():
    options, args = get_args()

    try:
        extract_vars(options)
    except IOError as err:
        sys.stderr.write("IOError " + str(err) + "\n")
        return


if __name__ == "__main__":
    sys.exit(main())

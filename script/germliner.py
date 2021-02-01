import sys, os
import fnmatch
from optparse import OptionParser
import pandas as pd
import ntpath
from collections import defaultdict
import json
import vcf
import time

pd.set_option('display.float_format', lambda x: '%.3f' % x)


def extract_all(options):
    dir = os.path.abspath(options.dir)
    print("Looking for .vcf files in directory: %s" % dir)

    for file in os.listdir(dir):
        if file.endswith("_snps_filt.vcf.gz"):
            options.freebayes_file = os.path.join(dir, file)
            print(options.freebayes_file)
            file_exists = os.path.join(os.getcwd(), ntpath.basename(options.freebayes_file).split("_")[0] + '_germline_snps.txt')

            if os.path.isfile(file_exists):
                print("Skipping sample %s" % ntpath.basename(options.freebayes_file).split("_")[0])
                continue

            time1 = time.time()
            sample, snps = parse_freebayes(options)
            time2 = time.time()
            write_vars(sample, snps, options)
            print('function took %0.3f s' % (time2-time1))

    return True


def find_normal(options):
    sample = ntpath.basename(options.freebayes_file).split("_")[0]
    config = pd.read_csv(options.config, delimiter="\t", index_col=False, na_filter=False, names=['sample', 'assay'])

    samples = config['sample'].tolist()
    it = iter(samples)

    s = {}
    for x in it:
        s[x] = next(it)

    if(s[sample]):
        print("Tumour: %s" % sample)
        print("Normal: %s" % s[sample])
        return sample, s[sample]
    else:
        print("Cant find corresponding normal sample for %s" % sample)


def write_vars(sample, snps, options):

    print("Total germline snps: %s" % (len(snps)))

    out_file = sample + '_germline_snps.txt'
    print("Writing snps to file %s" % out_file)

    header = '\t'.join(['chrom', 'pos', 'ref', 'alt', 'tumour_vaf', 'tumour_depth', 'normal_vaf', 'normal_depth'])

    with open(out_file, 'w') as snps_out:
        snps_out.write(header + '\n')

        for k in snps:
            chrom = snps[k]['info']['chrom']
            pos = snps[k]['info']['pos']
            alt = snps[k]['info']['alt']
            ref = snps[k]['info']['ref']

            taf = snps[k]['tumour']['VAF']
            tdp = snps[k]['tumour']['DP']

            naf = snps[k]['normal']['VAF']
            ndp = snps[k]['normal']['DP']

            l = '\t'.join(map(str, [chrom, pos, ref, alt, taf, tdp, naf, ndp]))
            snps_out.write(l + '\n')


def parse_freebayes(options):
    tumour, normal = find_normal(options)

    mode = 'r'
    if options.freebayes_file.endswith('.gz'):
        mode = 'rb'

    vcf_reader = vcf.Reader(open(options.freebayes_file, mode))
    chroms = ['2L', '2R', '3L', '3R', '4', 'X', 'Y']
    # chroms = 'Y'

    total_snp_count = 0
    snps = defaultdict(lambda: defaultdict(dict))

    # for c in chroms:
    #     print("Processing chrom: %s" % c)
    for record in vcf_reader:
        if record.CHROM not in chroms:
            continue

        if record.genotype(tumour)['DP'] and record.genotype(tumour)['DP'] > 20 and record.genotype(normal)['DP'] and record.genotype(normal)['DP'] > 20 and record.genotype(tumour)['GQ'] > 1:
            key = '_'.join(map(str, [record.CHROM, record.POS, str(record.ALT[0])]))

            if record.genotype(tumour)['GT'] != record.genotype(normal)['GT']:
                # print("Not equal genotyoes: %s vs %s" % (record.genotype(tumour)['GT'], record.genotype(normal)['GT']))
                continue

            total_snp_count += 1

            taf = round((record.genotype(tumour)['AO'] / (record.genotype(tumour)['AO'] + record.genotype(tumour)['RO'])),2)
            naf = round((record.genotype(normal)['AO'] / (record.genotype(normal)['AO'] + record.genotype(normal)['RO'])),2)
            tumour_gt, normal_gt = (record.genotype(tumour)['GT'], record.genotype(normal)['GT'])
            tumour_dp, normal_dp = (record.genotype(tumour)['DP'], record.genotype(normal)['DP'])

            t_dict = {'GT': tumour_gt, 'VAF': taf, 'DP': tumour_dp}
            n_dict = {'GT': normal_gt, 'VAF': naf, 'DP': normal_dp}
            info = {'chrom': record.CHROM, 'pos': record.POS, 'ref': record.REF,'alt': str(record.ALT[0]) }

            snps[key] = {
                'info': info,
                'tumour': t_dict,
                'normal': n_dict,
            }

    return tumour, snps


def main():
    parser = OptionParser()
    parser.add_option("-d", "--dir", dest="dir", help="Directory for batch processing")
    parser.add_option("-f", "--freebayes_file", dest="freebayes_file", help="Freebayes VCF file")
    parser.add_option("--config", dest="config", action="store", help="mapping for tumour/normal samples")

    parser.set_defaults(config='/Users/Nick_curie/Desktop/script_test/alleleFreqs/data/samples.tsv')

    options, args = parser.parse_args()

    if options.freebayes_file is None and options.dir is None:
        parser.print_help()
        print

    if options.dir:
        extract_all(options)

    else:
        try:
            time1 = time.time()
            sample, snps = parse_freebayes(options)
            time2 = time.time()
            write_vars(sample, snps, options)
            print('function took %0.3f s' % (time2-time1))

        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n")
            return


if __name__ == "__main__":
    sys.exit(main())

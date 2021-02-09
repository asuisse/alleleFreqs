import sys, os
import fnmatch
from optparse import OptionParser
import pandas as pd
import numpy as np
import ntpath
from collections import defaultdict
import json

pd.set_option('display.float_format', lambda x: '%.3f' % x)


def get_sample_names(options):
    sample_names = {}
    names = pd.read_csv(options.config, delimiter="\t", index_col=False, na_filter=False, names=['sample', 'sample_short', 'sample_converted', 'sex', 'assay'])
    sample_names = dict(zip(names['sample'], names['sample_converted']))

    return sample_names


def write_vars(snps, info, options):
    out_file = 'snps_combined.txt'
    print("Writing snps to file %s" % out_file)

    header = '\t'.join(['chrom', 'pos', 'ref', 'alt', 'depth', 'genotype', 'samples', 'sharedby'])

    with open(out_file, 'w') as snps_out:
        snps_out.write(header + '\n')

        for k in snps:
            chrom, pos, ref, alt = k.split('_')
            samples = snps[k]
            sharedby = len(samples)
            samples = ', '.join(samples)
            genotype = 'germline_recurrent'
            depth = 'na'
            if sharedby == 1:
                depth = info[k]['depth']
                if info[k]['genotype'] == 'germline':
                    genotype = 'germline_private'
                else:
                    genotype = info[k]['genotype']
            elif info[k]['genotype'] != 'germline':
                genotype = info[k]['genotype'] + '_recurrent'

            l = '\t'.join(map(str, [chrom, pos, ref, alt, depth, genotype, samples, sharedby]))
            snps_out.write(l + '\n')


def get_vars(options):
    dir = os.path.abspath(options.dir)
    print("Looking for files in directory: %s" % dir)

    sample_names = {}
    if options.config:
        sample_names = get_sample_names(options)

    overlaps = defaultdict(list)
    info = defaultdict(dict)

    excluded_samples = ["B241R41-2",  "A373R7", "A512R17"]

    # excluded_samples = []

    for file in os.listdir(dir):
        if file.endswith("_snps.txt"):
            sample = ntpath.basename(file).split("_")[0]
            if sample in excluded_samples:
                print("Skipping sample %s" % sample)
                continue

            file_path = os.path.join(dir, file)
            overlaps, info = extract_vars(sample, sample_names, overlaps, info, file_path)

    return overlaps, info
    # print(json.dumps(overlaps, indent=4, sort_keys=True))


def extract_vars(sample, sample_names, snps, info, f):
    df = pd.read_csv(f, delimiter="\t", index_col=False, na_filter=False)

    print("Old sample name: %s" % sample)

    if sample_names:
        sample = sample_names[sample]

    print("New sample name: %s" % sample)
    for i, row in df.iterrows():
        key = '_'.join(map(str, [row['chrom'], row['pos'], row['ref'], row['alt']]))
        snps[key].append(sample)
        info[key]['depth'] = np.mean([row['tumour_depth'], row['normal_depth']])
        info[key]['genotype'] = row['genotype']

    return snps, info


def main():
    parser = OptionParser()
    parser.add_option("-d", "--dir", dest="dir", help="Directory for batch processing")
    parser.add_option("-c", "--config", dest="config", help="Sample name mapping config")

    parser.set_defaults(dir=os.getcwd(), config='data/samples_names_conversion.txt')

    options, args = parser.parse_args()

    snps, info = get_vars(options)

    write_vars(snps, info, options)


if __name__ == "__main__":
    sys.exit(main())

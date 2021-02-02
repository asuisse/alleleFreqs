import sys, os
import fnmatch
from optparse import OptionParser
import pandas as pd
import ntpath
from collections import defaultdict
import json

pd.set_option('display.float_format', lambda x: '%.3f' % x)




def get_sample_names(options):
    sample_names = {}
    names = pd.read_csv(options.config, delimiter="\t", index_col=False, na_filter=False, names=['sample', 'sample_short', 'sample_converted', 'sex', 'assay'])
    sample_names = dict(zip(names['sample'], names['sample_converted']))

    return sample_names


def write_vars(snps, options):
    out_file = 'combined_snps.txt'
    print("Writing snps to file %s" % out_file)

    header = '\t'.join(['chrom', 'pos', 'ref', 'alt', 'samples', 'sharedby'])

    with open(out_file, 'w') as snps_out:
        snps_out.write(header + '\n')

        for k in snps:
            chrom, pos, ref, alt = k.split('_')
            samples = snps[k]
            sharedby = len(samples)
            samples = ', '.join(samples)
            status = 'recurrent'
            if sharedby == 1:
                status = 'private'

            l = '\t'.join(map(str, [chrom, pos, ref, alt, status, samples, sharedby]))
            snps_out.write(l + '\n')


def get_vars(options):
    dir = os.path.abspath(options.dir)
    print("Looking for files in directory: %s" % dir)

    sample_names = {}
    if options.config:
        sample_names = get_sample_names(options)

    overlaps = defaultdict(list)

    # excluded_samples = ["B241R41-2",  "A373R7", "A512R17", "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9", "D197R09", "D197R11", "D197R13", "D197R15", "D265R01", "D265R03", "D265R05", "D265R07", "D265R09", "D265R11", "D265R13"]

    excluded_samples = []

    for file in os.listdir(dir):
        if file.endswith("_germline_snps.txt"):
            sample = ntpath.basename(file).split("_")[0]
            if sample in excluded_samples:
                print("Skipping sample %s" % sample)
                continue

            file_path = os.path.join(dir, file)
            overlaps = extract_vars(sample, sample_names, overlaps, file_path)

    return overlaps
    # print(json.dumps(overlaps, indent=4, sort_keys=True))


def extract_vars(sample, sample_names, snps, f):
    df = pd.read_csv(f, delimiter="\t", index_col=False, na_filter=False)

    print("Old sample name: %s" % sample)

    if sample_names:
        sample = sample_names[sample]

    print("New sample name: %s" % sample)
    for i, row in df.iterrows():
        key = '_'.join(map(str, [row['chrom'], row['pos'], row['ref'], row['alt']]))
        snps[key].append(sample)

    return snps


def main():
    parser = OptionParser()
    parser.add_option("-d", "--dir", dest="dir", help="Directory for batch processing")
    parser.add_option("-c", "--config", dest="config", help="Sample name mapping config")

    parser.set_defaults(dir=os.getcwd(), config='data/samples_names_conversion.txt')

    options, args = parser.parse_args()

    snps = get_vars(options)

    write_vars(snps, options)


if __name__ == "__main__":
    sys.exit(main())

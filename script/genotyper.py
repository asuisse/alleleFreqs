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


def write_vars(snps, options):
    out_file = 'combined_snps.txt'
    print("Writing snps to file %s" % out_file)

    header = '\t'.join(['chrom', 'pos', 'ref', 'alt', 'samples'])

    with open(out_file, 'w') as snps_out:
        snps_out.write(header + '\n')

        for k in snps:
            chrom, pos, ref, alt = k.split('_')
            samples = snps[k]
            samples = ', '.join(samples)

            l = '\t'.join(map(str, [chrom, pos, ref, alt, samples]))
            snps_out.write(l + '\n')


def get_vars(options):
    dir = os.path.abspath(options.dir)
    print("Looking for files in directory: %s" % dir)

    overlaps = defaultdict(list)

    excluded_samples = ["B241R41-2",  "A373R7", "A512R17", "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9", "D197R09", "D197R11", "D197R13", "D197R15", "D265R01", "D265R03", "D265R05", "D265R07", "D265R09", "D265R11", "D265R13"]

    for file in os.listdir(dir):
        if file.endswith("_germline_snps.txt"):
            sample = ntpath.basename(file).split("_")[0]
            print(sample)
            if sample in excluded_samples:
                print("Skipping sample %s" % sample)
                continue
            overlaps = extract_vars(sample, overlaps, file)

    return overlaps
    # print(json.dumps(overlaps, indent=4, sort_keys=True))


def extract_vars(sample, snps, f):
    df = pd.read_csv(f, delimiter="\t", index_col=False, na_filter=False)

    for i, row in df.iterrows():
        key = '_'.join(map(str, [row['chrom'], row['pos'], row['ref'], row['alt']]))
        snps[key].append(sample)

    return snps


def main():
    parser = OptionParser()
    parser.add_option("-d", "--dir", dest="dir", help="Directory for batch processing")

    parser.set_defaults(dir=os.getcwd())

    options, args = parser.parse_args()

    snps = get_vars(options)

    write_vars(snps, options)


if __name__ == "__main__":
    sys.exit(main())

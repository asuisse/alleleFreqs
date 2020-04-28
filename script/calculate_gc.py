from __future__ import division
from __future__ import print_function

import fnmatch
import json
import ntpath
import os, re
import sys
from collections import defaultdict
from optparse import OptionParser

import pandas as pd
from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.Alphabet import IUPAC
# from Bio.SeqUtils import GC

from statistics import mean

import matplotlib.pyplot as plt
# from dna_features_viewer import BiopythonTranslator
import numpy as np

from scipy.interpolate import make_interp_spline, BSpline


def extract_loci(region):
    try:
        chrom, start, end = re.split('[:\\-]', region)
        chrom = chrom.replace('chr', '')
        start = start.replace(',', '')
        end = end.replace(',', '')

        print(chrom, start, end)
    except ValueError:
        sys.exit("[!] Must provide a valid region [ chrom:start:end ].Exiting.")

    return chrom, int(start), int(end)


def parse_variants(options):
    variants = pd.read_csv(options.variants, delimiter="\t", index_col=False, na_filter=False)
    plot = False

    if options.locus:
        chrom, start, end = extract_loci(options.locus)
        calc_gc(options, chrom, start, end, True, 'locus', 'locus')
        return True

    if options.event:
        variants = variants.loc[breakpoints['event'] == int(options.event)]
        plot = True

    sample = ntpath.basename(options.variants).split("_")[0]

    for index, row in variants.iterrows():
        if row['chromosome1'] == row['chromosome2']:
            calc_gc(options, row['chromosome1'], row['bp1'], row['bp2'], plot, row['event'], row['type'])

    return True


def parse_bed(options):
    regions = pd.read_csv(options.bedfile, delimiter="\t", index_col=False, na_filter=False, header=None, usecols=[0,1,2], names=['chrom', 'start', 'end'])
    plot = False

    gccontent = []

    for index, row in regions.iterrows():

        calc_gc_region(options, row['chrom'], row['start'], row['end'])

    # print("%.2f%%" % mean(gccontent))

    return True


def calc_gc_region(options, chrom, start, end):
    # genome = SeqIO.parse(options.genome, "fasta")
    record_dict = SeqIO.to_dict(SeqIO.parse(options.genome, "fasta"))

    print(record_dict['X'].seq[start:end])

    # for record in genome:

        # if record.id != chrom: continue
        # seq = record.seq[start:end]
        #
        # pos.append(w_start+i)
        # # print("-"*padding, seq, get_GC(seq))
        # padding += 1



    return True


def calc_gc(options, chrom, start, stop, plot, event, sv_type):
    genome = SeqIO.parse(options.genome, "fasta")

    if plot:
        w_start = int(start - 1e4)
        w_stop = int(stop + 1e4)
    else:
        w_start = start
        w_stop = stop

    print("Calculating GC content for event %s [%s bps %s] for windows of size %s over sequence: %s:%s-%s" % (event, sv_type, stop-start, options.window, chrom, w_start, w_stop))

    for record in genome:
        if record.id != chrom: continue
        percs = []
        pos = []
        if record.id == chrom:
            padding = 0
            for i in range(0, len(record.seq[w_start:w_stop]), int(options.window/2)):
                seq = record.seq[i:i+options.window]
                percs.append(get_GC(seq))
                pos.append(w_start+i)
                # print("-"*padding, seq, get_GC(seq))
                padding += 1

    print(" * Average GC%%: %s " % (int(mean(percs))))

    if plot: plot_gc_content(pos, percs, chrom, start, stop)

    return mean(percs)


def plot_gc_content(x, y, chrom, start, stop):
    x = [e / 1e6 for e in x]

    position = np.array(x)
    gc = np.array(y)
    xnew = np.linspace(position.min(),position.max(), 200)

    spl = make_interp_spline(position, gc, k=3)
    pos_smooth = spl(xnew)

    plt.plot(xnew,pos_smooth)
    plt.ylabel('GC %')
    plt.xlabel('Genomic position')
    plt.title('%% GC over region %s:%s-%s in %s' % ('2L', start, stop, 'sample'))

    plt.fill_between(xnew, pos_smooth, np.min(y), color='#539ecd')


    plt.ylim(bottom=0)
    plt.axvline(x=start/1e6, c='k')
    plt.axvline(x=stop/1e6, c='k')

    plt.show()


def get_GC(seq):
    c = 0
    for nucl in seq:
        if nucl in ['G', 'C']:
            c += 1
    return int(round(c/len(seq) * 100))


def get_args():
    parser = OptionParser()

    parser.add_option("-v", "--variants", dest="variants", action="store", help="svParser format file")
    parser.add_option("-b", "--bedfile", dest="bedfile", action="store", help="Bedfile containing regions")
    parser.add_option("-g", "--genome", dest="genome", action="store", help="Genome fasta file")
    parser.add_option("-w", "--window", dest="window", action="store", type="int", help="Window to slide over sequence. [Default 5]")
    parser.add_option("-l", "--locus", dest="locus", action="store", help="chr:start-end")
    parser.add_option("--event", dest="event", action="store", help="Inspect one event")

    parser.set_defaults(genome='/Users/Nick_curie/Documents/Curie/Data/Genomes/dmel_6.12.fa',
                        window=5)

    options, args = parser.parse_args()

    if not options.variants and not options.bedfile:
        parser.print_help()
        print()
        sys.exit("[!] Must provide a variants or bed file. Exiting.")
    else:
        return options, args


def main():
    options, args = get_args()
    try:
        if options.variants:
            parse_variants(options)
        if options.bedfile:
            parse_bed(options)
    except IOError as err:
        sys.stderr.write("IOError " + str(err) + "\n")
        return


if __name__ == "__main__":
    sys.exit(main())
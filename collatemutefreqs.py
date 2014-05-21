#!/usr/bin/env python

""" Utility to read the mutation frequency csvs from connor. """

import sys
import csv
import math
from scipy.stats import beta
from opener import opener
import utils

# TODO read versions from reference file to check

def fraction_uncertainty(n, total):
    """ Return uncertainty on the ratio n / total """
    assert n <= total
    n = float(n)  # in case they're integers
    total = float(total)
    if n == 0.0 or total == 0.0:
        return 0.0
    delta_n = math.sqrt(n)  # NOTE this isn't really right. For instance 1/2 will give negative lower bound
    delta_total = math.sqrt(total)
    delta_frac_squared = (n / total)**2 * ((delta_n/n)**2 + (delta_total/total)**2)
    return math.sqrt(delta_frac_squared)

class MuteFreqReader(object):
    def __init__(self, infname, imax=-1):
        self.maturities = ['memory', 'naive']  # NOTE eveywhere else I call this 'naivety' and give it the values 'M' or 'N'
        self.freqs = {}
        for human in utils.humans:
            self.freqs[human] = {}
        with opener('r')(infname) as infile:
            reader = csv.DictReader(infile)
            il = 0
            for line in reader:
                il += 1
                human = line['subject']
                gene_name = line['reference']
                if gene_name not in self.freqs[human]:
                    self.freqs[human][gene_name] = {}
                maturity = line['subset']
                if maturity not in self.freqs[human][gene_name]:
                    self.freqs[human][gene_name][maturity] = {}
                position = int(line['position'])
                # if position in self.freqs[human][gene_name][maturity]:
                #     print '  %d already in %s %s %s' % (position, human, gene_name, maturity)
                #     assert False
                assert position not in self.freqs[human][gene_name][maturity]
                self.freqs[human][gene_name][maturity][position] = {}
                self.freqs[human][gene_name][maturity][position]['ref'] = line['ref_base']
                self.freqs[human][gene_name][maturity][position]['n_reads'] = int(line['n_reads'])
                for nuke in utils.nukes:
                    self.freqs[human][gene_name][maturity][position][nuke] = float(line[nuke]) / int(line['n_reads'])
                if imax > 0 and il > imax:
                    break
    
    def print_freqs(self, only_human='', only_gene_name='', only_maturity='', only_position=-1):
        for human in self.freqs:
            if only_human != '' and human != only_human:
                continue
            print '%s -------' % human
            for gene_name in self.freqs[human]:
                if only_gene_name != '' and gene_name != only_gene_name:
                    continue
                for maturity in self.maturities:
                    if only_maturity != '' and maturity != only_maturity:
                        continue
                    print '  %-20s %-10s' % (gene_name, maturity)
                    mute_freqs = self.freqs[human][gene_name][maturity]
                    for position in sorted(mute_freqs):
                        if only_position > 0 and position != only_position:
                            continue
                        print_str = '    %3d %3s ' % (position, mute_freqs[position]['ref'])
                        for nuke in utils.nukes:
                            n = mute_freqs[position][nuke] * mute_freqs[position]['n_reads']
                            assert mute_freqs[position]['n_reads'] > 0
                            assert n >= 0
                            uncert = fraction_uncertainty(n, mute_freqs[position]['n_reads'])
                            print_str += '  | %4d / %4d = %8.3f +/- %3.3f' % (n, mute_freqs[position]['n_reads'], mute_freqs[position][nuke], uncert)
                        print print_str
                    
    
    

if __name__ == "__main__":
    for pair in [(1,1), (1,2), (3,3), (5,19), (50,190)]:
        print '%d/%d = %f +/- %f' % (pair[0], pair[1], float(pair[0])/pair[1], fraction_uncertainty(pair[0], pair[1]))
    sys.exit()
    human = 'A'
    naivety = 'M'
    mfr = MuteFreqReader('/home/dralph/Dropbox/work/recombinator/data/human-beings/' + human + '/' + naivety + '/mute-counts.csv.bz2', imax=1000)
    mfr.print_freqs(only_human='A', only_gene_name='IGHV1-NL1*01', only_maturity='memory')
    # mfr.print_freqs(only_human='C', only_gene_name='IGHV7-4-1*91', only_maturity='memory')


#!/usr/bin/env python

""" Utility to read the mutation frequency csv from connor. """

import sys
import csv
from opener import opener
import utils

# TODO read versions from reference file to check

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
                for maturity in self.maturities:
                    if maturity not in self.freqs[human][gene_name]:
                        self.freqs[human][gene_name][maturity] = {}
                    position = int(line['position'])
                    if position in self.freqs[human][gene_name][maturity]:
                        print '  %d already in %s %s %s' % (position, human, gene_name, maturity)
                        assert False
                    assert position not in self.freqs[human][gene_name][maturity]
                    self.freqs[human][gene_name][maturity][position] = {}
                    assert line['ref_base_memory'] == line['ref_base_naive']
                    self.freqs[human][gene_name][maturity][position]['ref'] = line['ref_base_naive']
                    for nuke in utils.nukes:
                        self.freqs[human][gene_name][maturity][position][nuke] = float(line[nuke + '_' + maturity])
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
                        print '    %3d %3s  %8f  %8f  %8f  %8f' % (position,
                                                                 mute_freqs[position]['ref'],
                                                                 mute_freqs[position]['A'],
                                                                 mute_freqs[position]['C'],
                                                                 mute_freqs[position]['G'],
                                                                 mute_freqs[position]['T'])
                    
    
    

if __name__ == "__main__":
    mfr = MuteFreqReader('/shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/IGHJ_mut_by_site.csv')
    mfr.print_freqs(only_human='C', only_gene_name='IGHJ6*03_F', only_maturity='memory')

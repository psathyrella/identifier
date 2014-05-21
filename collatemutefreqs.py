#!/usr/bin/env python

""" Utility to read the mutation frequency csvs from connor. """

import sys
import os
import csv
import math
sys.argv.append('-b')
from ROOT import TH1F, TCanvas, kRed, gROOT, TLine
gROOT.Macro("plotting/root/MitStyleRemix.cc+")
from scipy.stats import beta
from opener import opener
import utils

# TODO read versions from reference file to check

def fraction_uncertainty(obs, total):
    """ Return uncertainty on the ratio n / total """
    assert obs <= total
    if total == 0.0:
        return 0.0
    lo = beta.ppf(1./6, 1 + obs, 1 + total - obs)
    hi = beta.ppf(1. - 1./6, 1 + obs, 1 + total - obs)
    if float(obs) / total < lo:  # if k/n very small (probably zero), take a one-sided c.i. with 2/3 the mass
        lo = 0.
        hi = beta.ppf(2./3, 1 + obs, 1 + total - obs)
    if float(obs) / total > hi:  # same deal if k/n very large (probably one)
        lo = beta.ppf(1./3, 1 + obs, 1 + total - obs)
        hi = 1.
    return (lo,hi)

class MuteFreqReader(object):
    def __init__(self, infname, imax=-1):
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
                # assert line['N'] == ''
                for nuke in utils.nukes:
                    self.freqs[human][gene_name][maturity][position][nuke] = float(line[nuke]) / int(line['n_reads'])
                if imax > 0 and il > imax:
                    break
    
    def write_freqs(self, baseoutdir, total_frequency=False, only_human='', only_gene_name='', only_maturity='', only_position=-1, ):
        cvn = TCanvas("cvn", "", 1700, 600)
        for human in self.freqs:
            if only_human != '' and human != only_human:
                continue
            print '%s -------' % human
            for gene_name in self.freqs[human]:
                if only_gene_name != '' and gene_name != only_gene_name:
                    continue
                for maturity in utils.maturities:
                    if only_maturity != '' and maturity != only_maturity:
                        continue
                    print '  %-20s %-10s' % (gene_name, maturity)
                    mute_freqs = self.freqs[human][gene_name][maturity]
                    sorted_positions = sorted(mute_freqs)
                    hist = TH1F(human + '-' + utils.sanitize_name(gene_name) + '-' + maturity, '',
                                sorted_positions[-1] - sorted_positions[0] + 1,
                                sorted_positions[0] - 0.5, sorted_positions[-1] + 0.5)
                    lo_err_hist = TH1F(hist)
                    hi_err_hist = TH1F(hist)
                    for position in sorted_positions:
                        if only_position > 0 and position != only_position:
                            continue
                        # print_str = '    %3d %3s ' % (position, mute_freqs[position]['ref'])
                        n_conserved, n_mutated = 0, 0
                        total = mute_freqs[position]['n_reads']
                        for nuke in utils.nukes:
                            obs = int(round(mute_freqs[position][nuke] * total))
                            if nuke == mute_freqs[position]['ref']:
                                n_conserved += obs
                            else:
                                n_mutated += obs
                            uncert = fraction_uncertainty(obs, total)  # NOTE this is kinda slow
                            # print_str += ' %8.2f - %3.2f + %3.2f' % (mute_freqs[position][nuke], uncert[0], uncert[1])
                        if n_mutated + n_conserved != total:
                            print n_mutated,n_conserved,total
                            assert False
                        mutated_fraction = float(n_mutated) / total
                        mutated_fraction_err = fraction_uncertainty(n_mutated, total)
                        hist.SetBinContent(hist.FindBin(position), mutated_fraction)
                        lo_err_hist.SetBinContent(hist.FindBin(position), mutated_fraction_err[0])
                        hi_err_hist.SetBinContent(hist.FindBin(position), mutated_fraction_err[1])
                        # outfile.write('%20.12f %20.12f %20.12f %20.12f\n' % (position, mutated_fraction, mutated_fraction_err[0], mutated_fraction_err[1]))
                        # outfile.write(print_str + '\n')
                    hframe = TH1F(hist)
                    hframe.SetTitle(gene_name + ';;')
                    hframe.Reset()
                    hframe.SetMinimum(lo_err_hist.GetMinimum() - 0.03)
                    hframe.SetMaximum(1.1*hi_err_hist.GetMaximum())
                    hframe.Draw('')
                    line = TLine(hist.GetXaxis().GetXmin(), 0., hist.GetXaxis().GetXmax(), 0.)
                    line.SetLineColor(0)
                    line.Draw()  # can't figure out how to convince hframe not to draw a horizontal line at y=0, so... cover it up
                    hist.SetLineColor(419)
                    hist.Draw('same')
                    lo_err_hist.SetLineColor(kRed+2)
                    hi_err_hist.SetLineColor(kRed+2)
                    lo_err_hist.SetMarkerColor(kRed+2)
                    hi_err_hist.SetMarkerColor(kRed+2)
                    lo_err_hist.SetMarkerStyle(22)
                    hi_err_hist.SetMarkerStyle(23)
                    lo_err_hist.Draw('p same')
                    hi_err_hist.Draw('p same')
                    outdir = baseoutdir + '/' + human + '/' + utils.maturity_to_naivety(maturity) + '/plots'
                    if not os.path.exists(outdir):
                        os.makedirs(outdir)
                    outfname = outdir + '/' + utils.sanitize_name(gene_name) + '.png'
                    cvn.SaveAs(outfname)

if __name__ == "__main__":
    human = 'B'
    naivety = 'M'
    mfr = MuteFreqReader('/home/dralph/Dropbox/work/recombinator/data/human-beings/' + human + '/' + naivety + '/mute-counts.csv.bz2')
    mfr.write_freqs(os.getenv('PWD').replace('/home/dralph/Dropbox', os.getenv('www')), only_maturity='memory')
    # mfr.write_freqs(os.getenv('PWD').replace('/home/dralph/Dropbox', os.getenv('www')), only_human='A', only_gene_name='IGHV1-NL1*01', only_maturity='memory')

#!/usr/bin/env python
import csv
import os
import sys

from opener import opener
import utils
from Searcher import Searcher

with opener('r')('head-simulated-seqs.csv') as infile:
    germlines = utils.read_germlines('../../../recombinator')
    reader = csv.DictReader(infile)
    for inline in reader:
        print 'searching'
        searcher = Searcher(inline, debug=True, n_matches_max=10)
        searcher.search()
        print 'FOOP  ',
        for region in utils.regions:
            print '(',
            for imatch in range(len(searcher.matches[region])):
                print '%4.1e' % searcher.matches[region][imatch]['evalue'],
                if imatch > 2:
                    break
            print ') ',
        print ''

#!/usr/bin/env python
import csv
import os
import sys

from opener import opener
import utils
from Searcher import Searcher

values = {}
variables = ['found_strings', 'evalue', 'ali_from', 'ali_to', 'vd_insertion_length', 'dj_insertion_length']
for var in variables:
    values[var] = {}
    for region in utils.regions:
        values[var][region] = {}
        for imatch in range(3):
            values[var][region][imatch] = []
            
assert len(sys.argv) == 3;
data_type = sys.argv[1]
human = sys.argv[2]
# for human in A B C; do
#     datadir=data/human-beings/$human/M/data
#     bzgrep -m100 . $datadir/data.tsv.bz2 | sed 's/[ \t][ \t]*/,/g'|cut -f2 -d, |sed 's/nucleotide/seq/'> $datadir/head-data.csv
# done
naivety = 'M'
infname = ''
if data_type == 'simu':
    infname = '/home/dralph/Dropbox/work/recombinator/output/' + human + '/' + naivety + '/simu.csv'
else:
    infname = 'data/human-beings/' + human + '/' + naivety + '/' + data_type + '/head-data.csv'
baseoutdir = 'data/human-beings/' + human + '/' + naivety + '/' + data_type

print 'opening ',infname
print '  output',baseoutdir
with opener('r')(infname) as infile:
    germlines = utils.read_germlines('../../../recombinator')
    reader = csv.DictReader(infile)
    il = 0
    for inline in reader:
        il += 1
        print inline['seq'][-100:]
       # if len(inline['seq']) != 130:
       #     assert 'simulated' in infname
        searcher = Searcher(inline['seq'][-100:], debug=False, n_matches_max=5)
        found_str = searcher.search()
        values['found_strings']['v'][0].append(found_str)  # toss them in ['v'][0] -- doesn't really make sense, but they're fine anywhere
        if found_str != 'vjd':  # skip the ones where we didn't find matches in this order (see freqs above).
            continue
        for region in utils.regions:
            for imatch in range(len(searcher.matches[region])):
                if imatch > 2:
                    break
                match = searcher.matches[region][imatch]
                if imatch == 0 and region == 'd':
                   # print '%s (%3d%3d) --> (%3d%3d %s)' % (region, match['ali_from'], match['ali_to'], match['ali_from'] - 1, len(searcher.query_seqs[region]) - match['ali_to'], searcher.query_seqs[region]),
                    values['vd_insertion_length'][region][imatch].append(match['ali_from'] - 1)  # NOTE these are index *one* counting (!!!)
                    values['dj_insertion_length'][region][imatch].append(len(searcher.query_seqs[region]) - match['ali_to'])  # NOTE these are index *one* counting (!!!)
                values['evalue'][region][imatch].append(match['evalue'])
                values['ali_from'][region][imatch].append(match['ali_from'])
                values['ali_to'][region][imatch].append(match['ali_to'])
        # if il > 100:
        #     sys.exit()
       # print ''
       # break

for region in utils.regions:
    for var in variables:
        outdir = baseoutdir + '/' + var
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        for imatch in range(len(searcher.matches[region])):
            if imatch > 2:
                break
            with opener('w')(outdir + '/' + region + '-' + str(imatch) + '.txt') as outfile:
                for value in values[var][region][imatch]:
                    if var == 'found_strings':
                        outfile.write('%s\n' % value)
                    else:
                        outfile.write('%7e\n' % value)

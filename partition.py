#!/usr/bin/env python
import csv
import os

from opener import opener
import partutils
from Searcher import Searcher

def build_hmms():
    stodir = '/home/dralph/work/recombinator/data/msa'
    outdir = hmmdir
    for region in regions:
        for sto_fname in glob(stodir + '/' + region + '/*.sto'):
            if not os.path.exists(outdir + '/' + region):
                os.makedirs(outdir + '/' + region)
#            call(["../binaries/hmmbuild", '--dna', outdir + '/' + region + '/' + os.path.basename(sto_fname).replace('.sto','.hmm'), sto_fname]) #, stdout=foo)

    all_lines = []
    for region in regions:
        region_lines = []
        for hmmfname in glob(outdir + '/' + region + '/*.hmm'):
            with opener('r')(hmmfname) as hmmfile:
                region_lines.extend(hmmfile.readlines())
        with opener('w')(outdir + '/' + region + '/all.hmm') as pressfile:
            pressfile.writelines(region_lines)
        call(['../binaries/hmmpress', '-f', outdir + '/' + region + '/all.hmm'])
        all_lines.extend(region_lines)

    with opener('w')(outdir + '/all.hmm') as pressfile:
        pressfile.writelines(all_lines)
    call(['../binaries/hmmpress', '-f', outdir + '/all.hmm'])


# d gene freqs:
#    379 IGHD4-17*01
#    310 IGHD4-23*01
#     60 IGHD3-16*02
#     50 IGHD3-16*01
#     40 IGHD7-27*01
#     40 IGHD3-10*01
#     40 IGHD2/OR15-2b*01
#     40 IGHD2/OR15-2a*01
#     20 IGHD1-20*01
#     10 IGHD3-9*01
#     10 IGHD1-1*01

default_d = partutils.sanitize_name('IGHD4-17*01')  # it's the most common one, so use it if you can't figure out which it was

#build_hmms()
with opener('r')('head-simulated-seqs.csv') as infile:
    reader = csv.DictReader(infile)
    for inline in reader:
        searcher = Searcher(inline)
        searcher.search()
        inferred_group_str = ''
        true_group_str = ''
        for region in partutils.regions:
            inferred_name = searcher.get_best_match_name(region)
            true_name = partutils.sanitize_name(inline[region + '_gene'])
#            if inferred_name = 'none':
#                inferred_group_str += default_d
#                else:

            inferred_group_str += inferred_name
            true_group_str += true_name
            if inferred_name == 'none':
                print ' none',
            elif  inferred_name == true_name:
                print '  -  ',
            else:
                print '  x  ',
        for region in partutils.regions:
            print '%3d' % searcher.n_tries[region],
        print ''

        with opener('a')('inferred.csv') as outfile:
            print '  inferred: %s %d' % (inferred_group_str, hash(inferred_group_str))
            outfile.write(str(hash(inline['seq'])) + ',' + str(hash(inferred_group_str)) + '\n')
        with opener('a')('true.csv') as outfile:
            print '      true: %s %d' % (true_group_str, hash(true_group_str))
            outfile.write(str(hash(inline['seq'])) + ',' + str(hash(true_group_str)) + '\n')

# eg if [0] is 1e-5 and [1] is 1e-4, ratio is < 1
# or, if first guess is correct, its evalue should be very small, so ratio is very small
# whereas if they have similar evalues, ratio would be near one
# if ratio is greater than 1, the list wasn't sorted

# USE A BRT!

# So: if ratio > threshhold, think about using the second match
#   ... er, but what do I do with it? All I know is that it's *closer*
#   to being similar, not that it's *right*.
# turns out like this:
#     n_ok: 424 n_no: 75
#     ok: 1.21e-01 no: 2.89e-01
# so on the ones I get correct, the ratio is *smaller*, a factor two or so.

# oh, right, I guess what I was going to do was if the v id sucked (say, ratio greater than 1.75), identify the j first, chop
# it out, and see if that improves things. Or, even better, try each v that rates well against combinations of ds and js. We
# want the whole combo to rate well, after all.

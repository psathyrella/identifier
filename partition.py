#!/usr/bin/env python
import csv
import os
import sys

from opener import opener
import utils
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

#----------------------------------------------------------------------------------------
#build_hmms()
infname = 'head-simulated-seqs.csv'
#infname = 'head-data.csv'  # bzgrep -m100 . /shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/C/01-C-N_merged.tsv.bz2 | sed 's/[ \t][ \t]*/,/g'|cut -f2 -d, > head-data.csv
with opener('r')(infname) as infile:
    germlines = utils.read_germlines('../../../recombinator')
    reader = csv.DictReader(infile)
    for inline in reader:
        print 'searching'
#        inline['seq'] = inline['seq'][-130:]
        searcher = Searcher(inline['seq'], debug=True, n_matches_max=2)
        searcher.search()
        inferred_group_str = ''
        true_group_str = ''
        outline = {}
        outline['seq'] = inline['seq']
        print 'RESULT ',
        for region in utils.regions:
            inferred_name = searcher.get_best_match_name(region)
            outline[region + '_gene'] = utils.unsanitize_name(inferred_name)
            true_name = utils.sanitize_name(inline[region + '_gene'])

            inferred_group_str += inferred_name
            true_group_str += true_name
            if inferred_name == 'none':
                print ' none',
            elif  inferred_name == true_name:
                print '  -  ',
            else:
                print '  x  ',
        for region in utils.regions:
            print '%3d' % searcher.n_tries[region],
        print ''
        print '  true'
        utils.print_reco_event(germlines, inline, -1, -1)
        if searcher.all_matched():
            print '  inferred'
            try:
                searcher.build_inferred_seq(inline['seq'], germlines, outline)
                utils.print_reco_event(germlines, outline, -1, -1)
            except:
                print '   *something* is wrong!'
                print '    ',searcher.best_matches['v']
                print '    ',searcher.best_matches['d']
                print '    ',searcher.best_matches['j']
                continue
        else:
            print 'no matches!'
            print '    ',searcher.best_matches['v']
            print '    ',searcher.best_matches['d']
            print '    ',searcher.best_matches['j']
#            sys.exit()

        break

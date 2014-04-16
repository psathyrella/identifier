#!/usr/bin/env python
import csv
from subprocess import call,check_output
from glob import glob
import os
import sys
import math
import itertools

from opener import opener

regions = ('v','d','j')

hmmdir = 'data/hmms'

# set up a dict to access the columns
hmmer_output_column_list = ['total_evalue', 'total_score', 'total_bias', 'domain_evalue', 'domain_score', 'domain_bias', 'exp', 'N', 'Model']
hmmer_columns = {}
for icol in range(len(hmmer_output_column_list)):
    hmmer_columns[hmmer_output_column_list[icol]] = icol
    
#def fill_hmmer_output_line(line):
#    """ Takes a line of hmmer output, puts them in a list, and returns it. """
#    values = []
#    for column in hmmer_output_column_list:
#        values.append()
#
def sanitize_name(name):
    """ Replace characters in gene names that make crappy filenames. """
    saniname = name.replace('*', '_star_')
    saniname = saniname.replace('/', '_slash_')
    return saniname

def unsanitize_name(name):
    """ Re-replace characters in gene names that make crappy filenames. """
    unsaniname = name.replace('_star_', '*')
    unsaniname = unsaniname.replace('_slash_', '/')
    return unsaniname

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

def find_genes_in_seq(hmm_dbase_file, inline):
    seqfname = 'tmp/tmp.fa'
    with opener('w')(seqfname) as seqfile:
        seqlines = ['> tmp seq file\n', inline['seq'] + '\n']
        seqfile.writelines(seqlines)
    binary = "../binaries/hmmscan"
    # first cache the whole output in a file
    logfname = '/tmp/hmmer-log.txt'
#    call([binary, '--tblout ' + logfname, hmm_dbase_file, seq_file])  # wtf? the freaking *manual* says that's an option
    call([binary, '-o' + logfname, hmm_dbase_file, seqfname])
    # then filter it with sed to get the main block of scores
    filter_command = "sed -e '/>>/,/\[ok\]/ d' -e '/^#/ d' -e '/^$/ d' " + logfname + " | grep IGH"
    output = ''
    matches = []
    try:
        output = check_output(filter_command, shell=True)
    except:
        return matches
    for outline in output.splitlines():
        gene = unsanitize_name(outline.split()[hmmer_columns['Model']])
        score = float(outline.split()[hmmer_columns['total_score']])
        evalue = float(outline.split()[hmmer_columns['total_evalue']])
        matches.append({'gene': gene, 'score': score, 'evalue': evalue})

    return matches
    

#build_hmms()
#with opener('r')('/home/dralph/work/recombinator/out.csv') as infile:
with opener('r')('out.csv') as infile:
    reader = csv.DictReader(infile)
    for region in regions:
        print '%s region' % region
        total_ratio = {}
        total_ratio['ok'] = 0.0
        total_ratio['no'] = 0.0
        n_tot = {}
        n_tot['ok'] = 0
        n_tot['no'] = 0
        for inline in reader:
            matches = find_genes_in_seq(hmmdir + '/' + region + '/all.hmm', inline)
            simulated_gene = inline[region + '_gene']
            if len(matches) == 0:
                print '  no match'
                continue

            match_first = 'no'  # does the correct gene have the highest score (lowest evalue)?
            if simulated_gene == matches[0]['gene']:
                match_first = 'ok'
            print '  %s' % match_first,
            for im in range(1, len(matches)):  # look at the next few matches
                if im > 1:  # only look at the second match a.t.m.
                    break
                ratio = matches[0]['evalue'] / matches[im]['evalue']
                total_ratio[match_first] += ratio
                n_tot[match_first] += 1
                match_this = 'no'  # does the correct gene have the highest score (lowest evalue)?
                if simulated_gene == matches[im]['gene']:
                    match_this = 'ok'
                print ' (' + match_this + ',',
                print '%6.0e)' % ratio,
            print ''
        print '   n_ok: %d n_no: %d' % (n_tot['ok'], n_tot['no'])
        print '   ok: %6.2e no: %6.2e' % (total_ratio['ok'] / n_tot['ok'], total_ratio['no'] / n_tot['no'])
        sys.exit()

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

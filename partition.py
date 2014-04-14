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
hmmer_output_column_list = ['total_E-value', 'total_score', 'total_bias', 'domain_E-value', 'domain_score', 'domain_bias', 'exp', 'N', 'Model']
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

#build_hmms()
with opener('r')('/home/dralph/work/recombinator/out.csv') as infile:
    reader = csv.DictReader(infile)
    for inline in reader:
        tmpfile = '/tmp/tmp.fa'
        with opener('w')(tmpfile) as seqfile:
            seq_lines = ['> tmp seq file\n', inline['seq'] + '\n']
            seqfile.writelines(seq_lines)
#        output = check_output(["../binaries/hmmscan", hmmdir + '/all.hmm', tmpfile])
        binary = "../binaries/hmmscan"
        hmm_dbase = hmmdir + '/all.hmm'
        seq_file = tmpfile
        logfname = '/tmp/hmmer-log.txt'
        # first cache the whole output in a file
        with opener('w')(logfname) as logfile:
            call([binary, hmm_dbase, seq_file], stdout=logfile)
        # then filter it with sed to get the main block of scores
        filter_command = "sed -e '/>>/,/\[ok\]/ d' -e '/^#/ d' -e '/^$/ d' " + logfname + " | grep IGH"
        output = check_output(filter_command, shell=True)
        matches = {}
        for region in regions:
            matches[region] = []
        for outline in output.splitlines():
            gene_name = outline.split()[hmmer_columns['Model']]
            region = gene_name[3].lower()
            assert region=='v' or region =='d' or region=='j'
            matches[region].append(outline)
        warn_string = ''
        for region in regions:
            if len(matches[region]) == 0:
                print 'na',
            for match in matches[region]:
                # make sure almost all the match score is from one domain
                #print match.split()[hmmer_columns['total_score']]
                domain_score = float(match.split()[hmmer_columns['domain_score']])
                total_score = float(match.split()[hmmer_columns['total_score']])
                simulated_gene = inline[region + '_gene']
                inferred_gene = unsanitize_name(match.split()[hmmer_columns['Model']])
                is_correct = simulated_gene == inferred_gene
                if is_correct:
                    print ' -',
                else:
                    print ' x',
                if domain_score / total_score < 0.9:
                    warn_string += '    WARNING: domain / total = %.2f / %.2f = %.2f for %s (%s)' % (domain_score, total_score, domain_score/total_score, simulated_gene, inferred_gene)
#                    # print more info about that sequence:
#                    with opener('r')(logfname) as logfile:
#                        for line in itertools.dropwhile(lambda line: '>> ' + gene_name not in line, logfile.readlines()):
#                        #for line in itertools.islice(logfile.readlines(), lambda line: '>> ' + gene_name not in line,):
#                            print line.strip('\n')
                break  # break after the best match (they should be sorted)
        print warn_string

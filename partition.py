#!/usr/bin/env python
import csv
from subprocess import call, check_call, check_output, Popen, PIPE, STDOUT
from itertools import ifilter, dropwhile, takewhile, tee
from glob import glob
import os
import sys
import math

from opener import opener

Hmmerer(region, '/home/dralph/Dropbox/work/hmmer/hmmer-3.1b1-linux-intel-x86_64')

regions = ('v','d','j')

hmmdir = 'data/hmms'

# set up a dict to access the columns
hmmer_output_column_list = ['target_name', 'accession', 'query_name', 'accession', 'hmm_from', 'hmm_to', 'ali_from', 'ali_to', 'env_from', 'env_to', 'mod_len', 'strand', 'evalue', 'score', 'bias', 'description']
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

def get_match_info(length, command_output, gene):
    """ Parse hmmer output for info associated with <gene>. """
    if length == 'short':
        print command_output
        sys.exit()
        matchlines = list(ifilter(lambda line: line.find(sanitize_name(gene)) >= 0, command_output.splitlines()))  # get the lines that have the gene name in them (there should only be one)
        assert len(matchlines) == 1
        return matchlines[0]
    elif length == 'long':
        # drop the pipeline summary crap
        matchlines = takewhile(lambda line: line.find('Internal pipeline statistics summary:') < 0, command_output.splitlines())
        # drop all lines until the '>> <gene>' line
        matchlines = dropwhile(lambda line: line.strip() != '>> ' + sanitize_name(gene), matchlines)
        # then drop the '>> <gene>' line
        matchlines = dropwhile(lambda line: line.strip() == '>> ' + sanitize_name(gene), matchlines)
        # then drop the lines after and including the '>> <next_gene>' line
        matchlines = takewhile(lambda line: line.find('>> IGH') < 0, matchlines)
        # get the summary info for each domain. first is the header line
        header_line = matchlines.next()
        if header_line.split()[3] != 'c-Evalue':  # make sure the column order hasn't changed
            print header_line,'---->',header_line.split()[3],'c-Evalue'  # make sure the column order hasn't changed
        assert header_line.split()[3] == 'c-Evalue'  # make sure the column order hasn't changed
        assert header_line.split()[8] == 'alifrom'
        assert header_line.split()[9] + header_line.split()[10] == 'alito'
        # TODO use c-Evalue or i-Evalue or a combo?
        matchlines.next()  # throw away the line of dashes
        domains = []
        for line in matchlines:
            if line.strip() == '':  # break when we get to the end of the domain summary section
                break
            vals = line.split()
            conditional_evalue = float(vals[4])
            independent_evalue = float(vals[5])
            ali_from = int(vals[9])
            ali_to = int(vals[10])
            domains.append({'conditional_evalue': conditional_evalue, 'independent_evalue':independent_evalue, 'ali_from': ali_from, 'ali_to': ali_to})

        return domains
    else:
        assert False

def find_matches(hmmdir, region, seq, n_tries, debug=False):
#    print ' %s' % region,
    if region not in n_tries:
        n_tries[region] = 0
    hmm_dbase_file = hmmdir + '/' + region + '/all.hmm'
    seqfname = '-'  # read query seq from stdin
    query_seq_command = 'echo \"> test seq file\n' + seq + '\n\"'  # pipe it to the binary to avoid writing to disk
    hmmerdir = '/home/dralph/Dropbox/work/hmmer/hmmer-3.1b1-linux-intel-x86_64'
    binary = '/bin/bash -c \"' + hmmerdir + '/binaries/nhmmscan'
    options = '--tblout >(cat)'
    if n_tries[region] > 0:
        options += ' --max'
#        binary = binary.replace('nhmmscan', 'hmmscan')  # docs are a little unclear on the real differences -- but in practice hmmscan gets more matches

    # run that sucker!
    output = ''
    n_tries[region] += 1
    # oh, the contortions! just to get decent output and not write it to disk
    command = query_seq_command + ' | ' + binary + ' ' + options + ' ' + hmm_dbase_file + ' ' + seqfname + '>/dev/null\"'
    if debug:
        print command
    output = check_output(command, shell=True)
    matches = []
    outlines = output.splitlines()
    header_line = outlines[0].split()
    # make sure column order hasn't changed
    assert header_line[1] + header_line[2] == 'targetname'
    assert header_line[7] == 'hmmfrom'
    assert header_line[8] + header_line[9] == 'hmmto'
    assert header_line[10] == 'alifrom'
    assert header_line[11] + header_line[12] == 'alito'
    assert header_line[18] == 'E-value'
    for outline in outlines[2:]:  # this would be a whole fucking lot easier if they would make their column headers one damn word each
        if outline[0] == '#':  # skip crap at the end
            break
        entry = {}
        entry['target_name'] = outline.split()[hmmer_columns['target_name']]
        for column in ['hmm_from', 'hmm_to', 'ali_from', 'ali_to']:
            entry[column] = int(outline.split()[hmmer_columns[column]])
        entry['evalue'] = float(outline.split()[hmmer_columns['evalue']])
        matches.append(entry)
        if len(matches) > 1:  # only look at the first two matches a.t.m.
            break

#    if len(matches) == 0:  # no matches
#        print ' (none)',
#    else:
#        print ' (' + matches[0]['target_name']+ ')',
    return matches
    
def excise_match(seq, matches):
    """ Once we're sure enough we have the correct germ line gene
    for one region {vdj}, cut it out to help with identification
    of the others.
    """

    i_best = -1
    best_match = {}
    if len(matches) == 1:
        i_best = 0
        best_match = matches[0]
    else:  # more than one match -- figure out what to do about that
        # start by finding the best match
        best_evalue = 999
        for imatch in range(len(matches)):
            match = matches[imatch]
            if match['evalue'] < best_evalue:
                best_evalue = match['evalue']
                i_best = imatch
                best_match = match

    # NOTE convert from 1-start to 0-start indexing
    excise_from = matches[i_best]['ali_from'] - 1
    excise_to = matches[i_best]['ali_to'] - 1
    seq = seq[:excise_from] + seq[excise_to + 1 :]
    return seq

def execute_search(hmmdir, region, seq, is_matched, n_tries, simulated_gene, is_correct, debug=False):
    """ Look for matches for <region>. If you find them, cut them out and return the shorter sequence. """
    matches= find_matches(hmmdir, region, seq, n_tries, debug=debug)
    if len(matches) != 0:
        is_matched[region] = True
        if sanitize_name(simulated_gene) == matches[0]['target_name']:
            is_correct[region] = '-'
        else:
            is_correct[region] = 'x'
        return excise_match(seq, matches)   # if there's a good match cut it out. TODO change 'good' criterion
    else:
        is_matched[region] = False
        is_correct[region] = 'none'
        return seq


#build_hmms()
with opener('r')('simulated-seqs.csv') as infile:
    reader = csv.DictReader(infile)
    for inline in reader:
        is_matched, n_tries, is_correct = {}, {}, {}
        current_seq = inline['seq']
        for region in ['v', 'j', 'd']:  # do d last
            current_seq = execute_search(hmmdir, region, current_seq, is_matched, n_tries, inline[region + '_gene'], is_correct)
#        print ''
        # now try again!
        for region in ['v', 'j', 'd']:  # do d last, still
            if not is_matched[region]:
                current_seq = execute_search(hmmdir, region, current_seq, is_matched, n_tries, inline[region + '_gene'], is_correct)  # NOTE current_seq potentially has bits cut out of it
            if not is_matched[region]:
                print '\n  dammit %s still not matched: %s' % (region, current_seq),
                execute_search(hmmdir, region, current_seq, is_matched, n_tries, inline[region + '_gene'], is_correct, debug=True)
        print '  %4s %4s %4s' % (is_correct['v'], is_correct['d'],is_correct['j'],)
#        print '\n----------------------------------------------------------------------------------------'
#        break

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

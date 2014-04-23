#!/usr/bin/env python
import csv
from subprocess import call, check_output, Popen, PIPE, STDOUT
from itertools import ifilter, dropwhile, takewhile, tee
from glob import glob
import os
import sys
import math

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

def get_match_info(length, command_output, gene):
    """ Parse hmmer output for info associated with <gene>. """
    if length == 'short':
        matchlines = list(ifilter(lambda line: line.find(sanitize_name(gene)) >= 0, command_output.splitlines()))  # get the lines that have the gene name in them (there should only be one)
        assert len(matchlines) == 1
        return matchlines[0]
    elif length == 'long':
#        print command_output
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
        # and skip the rest of the info for the moment
#        # and for the time being, drop any second domain information
#        # TODO does that really make sense
#        matchlines = takewhile(lambda line: line.find('== domain 2') < 0, matchlines)
#        # TODO make sure only one domain is listed
#        germline_seq, matching_position_line, test_seq, probabilities = '', '', '', ''
#        for line in matchlines:
#            # TODO I think this will fail if there's indels in the interior
#            if len(line) == 0:
#                continue
#            print line
#            # first get the germline line
#            line = line.split()
#            assert line[0] == sanitize_name(gene)
#            start_position = line[1]  # start position of this match
#            germline_seq += line[2]
#            end_position = line[3]
#            # then get the line with positions at which they match
#            matching_position_line +=  matchlines.next()  # screwed up with extra spaces at the start, but I'm not using it a.t.m.
#            # then the line for the test sequence
#            line = matchlines.next().split()
#            assert line[0] == 'test'
##            assert line[1] == start_position  # um, maybe not?
#            test_seq += line[2]
##            assert line[3] == end_position  # um, maybe not?
#            # and finally, get the probability line
#            probabilities += matchlines.next().split()[0]
#        return {'germline_seq': germline_seq, 'test_seq': test_seq}
    else:
        assert False

def find_matches(hmmdir, region, seq, n_tries):
    print ' %s' % region,
    if region not in n_tries:
        n_tries[region] = 0
    hmm_dbase_file = hmmdir + '/' + region + '/all.hmm'
    seqfname = 'tmp/tmp.fa'
    with opener('w')(seqfname) as seqfile:
        seqlines = ['> test seq file\n', seq + '\n']
        seqfile.writelines(seqlines)
    hmmerdir = '/home/dralph/Dropbox/work/hmmer/hmmer-3.1b1-linux-intel-x86_64'
    binary = hmmerdir + "/binaries/hmmscan"  # TODO read from stdin
    options = []
#    if region == 'd':  # for d region we really gotta scrape the bottom of the barrel
#        pass
#        options += ['--incE', '1000']
#        options += ['--incT', '0']
#        options += ['--T', '0']
    if n_tries[region] > 0:
        options += ['--max', ]
    full_output = ''
    try:
        full_output = check_output([binary,] + options + [hmm_dbase_file, seqfname], stderr=STDOUT)
        # TODO arg, trying to get the output if it fails. sigh
    except:
        print full_output
    # then filter it with sed to get the main block of scores
    filter_command = "sed -e '/>>/,/\[ok\]/ d' -e '/^#/ d' -e '/^$/ d' | grep IGH"  # TODO rewrite this with itertools
    proc = Popen(filter_command, shell=True, stdout=PIPE, stdin=PIPE)
    summary_output = (proc.communicate(input=full_output)[0]).strip()  # the block of output at the top with a line per match
    if len(summary_output) == 0:  # no matches
        print '  no matches'
        return []
    matches = []
    for outline in summary_output.splitlines():
        gene = unsanitize_name(outline.split()[hmmer_columns['Model']])
        score = float(outline.split()[hmmer_columns['total_score']])
        evalue = float(outline.split()[hmmer_columns['total_evalue']])
        matches.append({'gene': gene,
                        'score': score,
                        'evalue': evalue,
                        'short_info': get_match_info('short', summary_output, gene),
                        'long_info': get_match_info('long', full_output, gene)})
        if len(matches) > 1:  # only look at the first two matches a.t.m.
            break

    n_tries[region] += 1
    return matches
    
def excise_match(seq, gene_name, hmmer_output):
    """ Once we're sure enough we have the correct germ line gene
    for one region {vdj}, cut it out to help with identification
    of the others.
    """

    i_best = -1
    best_domain = {}
    if len(hmmer_output) == 1:
        i_best = 0
        best_domain = hmmer_output[0]
    else:  # more than one domain -- figure out what to do about that
        # start by finding the best domain
        best_evalue = 999
        for idom in range(len(hmmer_output)):
            domain = hmmer_output[idom]
            if domain['conditional_evalue'] < best_evalue:
                best_evalue = domain['conditional_evalue']
                i_best = idom
                best_domain = domain
        # now loop through the others to figure out what to do about them
        # aw, screw it, this is turning out kind of complicated to decide what to do about the other matches. For now just ignore them. TODO use the fact that you know the right end of j and left end of v won't be eroded
#        for idom in range(len(hmmer_output)):
#            if idom == i_best:
#                continue
#            domain = hmmer_output[idom]
#            ali_from = domain['ali_from']
#            ali_to = domain['ali_to']
#            # if they have non-overlapping regions, smash em together
#            # NOTE the evalues in best_domain will no longer be correct after this
#            if ali_to > best_domain['ali_to']:  # this domain extends beyond the best one
#                best_domain['ali_to'] = ali_to  # smash 'em together. NOTE that this will *include* any non-matched regions between them. which is what I want for vdj recombination
#            if ali_from < best_domain['ali_from']:  # this domain starts before the best one
#                best_domain['ali_from'] = ali_from  # smash 'em together. NOTE that this will *include* any non-matched regions between them. which is what I want for vdj recombination
#            # and if they don't... um... just ignore them for now. TODO is that the best way to handle this?

    # NOTE convert from 1-start to 0-start
    excise_from = hmmer_output[i_best]['ali_from'] - 1
    excise_to = hmmer_output[i_best]['ali_to'] - 1
#    print '    excising from %d to %d' % (excise_from, excise_to)
#    print '   ',seq
    seq = seq[:excise_from] + seq[excise_to + 1 :]
#    print '   ',seq
    return seq

def execute_search(hmmdir, region, seq, is_matched, n_tries):
    """ Look for matches for <region>. If you find them, cut them out and return the shorter sequence. """
    matches = find_matches(hmmdir, region, seq, n_tries)
    if len(matches) != 0:  # if there's a good v match cut it out. TODO change 'good' criterion
        is_matched[region] = True
        return excise_match(seq, matches[0]['gene'], matches[0]['long_info'])
    else:
        is_matched[region] = False
        return seq

#build_hmms()
with opener('r')('head-simulated-seqs.csv') as infile:
    reader = csv.DictReader(infile)
    for inline in reader:
        is_matched, n_tries = {}, {}
        current_seq = inline['seq']
        for region in ['v', 'j', 'd']:  # do d last
            current_seq = execute_search(hmmdir, region, current_seq, is_matched, n_tries)
        # now try again!
        for region in ['v', 'j', 'd']:  # do d last, still
            if not is_matched[region]:
                current_seq = execute_search(hmmdir, region, current_seq, is_matched, n_tries)  # NOTE current_seq potentially has bits cut out of it
            if not is_matched[region]:
                print '  dammit %s still not matched: %s' % (region, current_seq)
        print ''
#
#
#        
#        matches = {}
#        region = 'v'
#        matches[region] = find_matches(hmmdir, region, inline['seq'])
#        excised_seq = inline['seq']
#        if len(matches[region]) != 0:  # if there's a good v match cut it out. TODO change 'good' criterion
##            print '%s match %18s %8r' % (region, inline[region + '_gene'], inline[region + '_gene'] == matches[region][0]['gene'])
#            excised_seq = excise_match(inline['seq'], matches[region][0]['gene'], matches[region][0]['long_info'])
#        # so then try to find a j
#        region = 'j'
#        matches[region] = find_matches(hmmdir, region, excised_seq)
#        if len(matches[region]) != 0:
##            print 'j match %18s %8r' % (inline[region + '_gene'], inline[region + '_gene'] == matches[region][0]['gene'])
#            excised_seq = excise_match(excised_seq, matches[region][0]['gene'], matches[region][0]['long_info'])
#        # now, try to get a d!
#        region = 'd'
#        matches[region] = find_matches(hmmdir, region, excised_seq)
#        if len(matches[region]) != 0:
#            print 'd match %18s %8r' % (inline[region + '_gene'], inline[region + '_gene'] == matches[region][0]['gene'])

#        simulated_gene = inline[region + '_gene']
#        match_first = 'no'  # does the correct gene have the highest score (lowest evalue)?
#        if simulated_gene == matches[0]['gene']:
#            match_first = 'ok'
#        print '  %s' % match_first,
#        sys.exit()
#        for im in range(1, len(matches)):  # look at the next few matches
#            if im > 1:  # only look at the second match a.t.m.
#                break
##            ratio = matches[0]['evalue'] / matches[im]['evalue']
#        print ''

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

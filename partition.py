#!/usr/bin/env python
import csv
import os
import sys

from Bio import SeqIO

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
#default_d = partutils.sanitize_name('IGHD4-17*01')  # it's the most common one, so use it if you can't figure out which it was

def read_germlines():
    germlines = {}
    for region in partutils.regions:
        germlines[region] = {}
        for seq_record in SeqIO.parse('../../../recombinator/data/igh'+region+'.fasta', "fasta"):
            germlines[region][seq_record.name] = str(seq_record.seq)
    return germlines

def print_debug_info(seq, true_name, inferred_name, germlines):
    pass

class Colors:
    head = '\033[95m'
    blue = '\033[94m'
    green = '\033[92m'
    yellow = '\033[93m'
    red = '\033[91m'
    end = '\033[0m'

def red(seq):
    return Colors.red + seq + Colors.end

def is_mutated(original, final):
    if original == final:
        return final
    else:
        return red(final)

def print_reco_event(germlines, line):
    # NOTE took this code from Recombinator -- should break it out into a util fcn
    original_seqs = {}
    for region in partutils.regions:
        original_seqs[region] = germlines[region][line[region+'_gene']]
    v_start = 0
    v_length = len(original_seqs['v']) - int(line['v_3p_del'])
    d_start = v_length + len(line['vd_insertion'])
    d_length = len(original_seqs['d']) - int(line['d_5p_del']) - int(line['d_3p_del'])
    j_start = v_length + len(line['vd_insertion']) + d_length + len(line['dj_insertion'])
    j_length = len(original_seqs['j']) - int(line['j_5p_del'])
    eroded_seqs = {}
    eroded_seqs['v'] = original_seqs['v'][:len(original_seqs['v'])-int(line['v_3p_del'])]
    eroded_seqs['d'] = original_seqs['d'][int(line['d_5p_del']) : len(original_seqs['d'])-int(line['d_3p_del'])]
    eroded_seqs['j'] = original_seqs['j'][int(line['j_5p_del']) :]

    germline_v_end = len(original_seqs['v']) - 1
    germline_d_start = len(original_seqs['v']) - int(line['v_3p_del']) + len(line['vd_insertion']) - int(line['d_5p_del'])
    germline_d_end = germline_d_start + len(original_seqs['d'])
    germline_j_start = germline_d_end + 1 - int(line['d_3p_del']) + len(line['dj_insertion']) - int(line['j_5p_del'])

    final_seq = ''
    for inuke in range(len(line['seq'])):
        ilocal = inuke
        if ilocal < v_length:
            final_seq += is_mutated(eroded_seqs['v'][ilocal], line['seq'][inuke])
        else:
            ilocal -= v_length
            if ilocal < len(line['vd_insertion']):
                final_seq += is_mutated(line['vd_insertion'][ilocal], line['seq'][inuke])
            else:
                ilocal -= len(line['vd_insertion'])
                if ilocal < d_length:
                    final_seq += is_mutated(eroded_seqs['d'][ilocal], line['seq'][inuke])
                else:
                    ilocal -= d_length
                    if ilocal < len(line['dj_insertion']):
                        final_seq += is_mutated(line['dj_insertion'][ilocal], line['seq'][inuke])
                    else:
                        ilocal -= len(line['dj_insertion'])
                        final_seq += is_mutated(eroded_seqs['j'][ilocal], line['seq'][inuke])

    # pad with dots
    eroded_seqs['v'] = eroded_seqs['v'] + int(line['v_3p_del']) * '.'
    eroded_seqs['d'] = int(line['d_5p_del']) * '.' + eroded_seqs['d'] + int(line['d_3p_del']) * '.'
    eroded_seqs['j'] = int(line['j_5p_del']) * '.' + eroded_seqs['j']

    insertions = v_length * ' ' + line['vd_insertion'] + d_length * ' ' + line['dj_insertion'] + j_length * ' '
    vj = germline_d_start * ' ' + eroded_seqs['d'] + (len(original_seqs['j']) - int(line['j_5p_del']) + len(line['dj_insertion']) - int(line['d_3p_del'])) * ' '
    d = eroded_seqs['v'] + (germline_j_start - germline_v_end - 2) * ' ' + eroded_seqs['j']

    print '    ',insertions
    print '    ',final_seq
    print '    ',vj
    print '    ',d
    if 'ack' in line and line['ack']:
        sys.exit()
#    assert len(line['seq']) == line['v_5p_del'] + len(hmms['v']) + len(outline['vd_insertion']) + len(hmms['d']) + len(outline['dj_insertion']) + len(hmms['j']) + outline['j_3p_del']

#build_hmms()
with opener('r')('head-simulated-seqs.csv') as infile:
    germlines = read_germlines()
    reader = csv.DictReader(infile)
    for inline in reader:
        print 'searching'
        searcher = Searcher(inline, debug=True)
        searcher.search()
        inferred_group_str = ''
        true_group_str = ''
        outline = {}
        outline['seq'] = inline['seq']
        for region in partutils.regions:
            inferred_name = searcher.get_best_match_name(region)
            outline[region + '_gene'] = partutils.unsanitize_name(inferred_name)
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
        print '  true'
        print_reco_event(germlines, inline)
        searcher.build_inferred_seq(inline['seq'], germlines, outline)
#        build_line_from_inferred_names(inline['seq'], outline, germlines, searcher)
        print '  inferred'
        print_reco_event(germlines, outline)
#        for key,val in outline.iteritems():
#            print key,val
#        print_debug_info(true_name, inferred_name, germlines)

#        break
#        with opener('a')('inferred.csv') as outfile:
#            print '  inferred: %s %d' % (inferred_group_str, hash(inferred_group_str))
#            outfile.write(str(hash(inline['seq'])) + ',' + str(hash(inferred_group_str)) + '\n')
#        with opener('a')('true.csv') as outfile:
#            print '      true: %s %d' % (true_group_str, hash(true_group_str))
#            outfile.write(str(hash(inline['seq'])) + ',' + str(hash(true_group_str)) + '\n')

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

""" Execute search with proper hmmer binary, and store output in a dict. """

import sys
from subprocess import check_call,check_output
from itertools import dropwhile

#----------------------------------------------------------------------------------------
def sanitize_header_line(line):
    """ Really, dude? You couldn't make the damn headers one word each? """
    line = line.lstrip('#')
    line = line.replace('target name', 'target_name')
    line = line.replace('query name', 'query_name')
    line = line.replace('hmmfrom', 'hmm_from')
    line = line.replace('hmm to', 'hmm_to')
    line = line.replace('alifrom', 'ali_from')
    line = line.replace('ali to', 'ali_to')
    line = line.replace('envfrom', 'env_from')
    line = line.replace('env to', 'env_to')
    line = line.replace('E-value', 'evalue')
    line = line.replace('description of target', 'description')
    return line

#----------------------------------------------------------------------------------------
def sanitize_hmmscan_header_line(line):
    """ *and* you had to make it a *two* line header for hmmscan? """
    # NOTE this assumes the order is hmm, ali, env
    line = line.replace(' from', 'hmm_from', 1)
    line = line.replace(' to', 'hmm_to', 1)
    line = line.replace(' from', 'ali_from', 1)
    line = line.replace(' to', 'ali_to', 1)
    line = line.replace(' from', 'env_from', 1)
    line = line.replace(' to', 'env_to', 1)
    return line

#----------------------------------------------------------------------------------------
class Hmmerer(object):
    """ ya what it says up there yo """

    def __init__(self, hmmerdir, region, seq, sensitivity=''):
        self.region = region
        self.seq = seq
        self.sensitivity = sensitivity
        self.binary = 'nhmmscan'
        self.command = self.build_command(hmmerdir)
        self.output = ''  # raw output from hmmer
        self.matches = []

    def build_command(self, hmmerdir):
        """ Build hmmer command. """
        hmm_dbase_file = 'data/hmms/' + self.region + '/all.hmm'
        query_seq_command = 'echo \"> test seq file\n' + self.seq + '\n\"'  # command to generate query sequence 'file': pipe it to hmmer to avoid writing to disk
        options = ''
        if self.sensitivity == 'more':
            options += ' --max -E 1000 --incE 1000 '
        elif self.sensitivity == 'max':
            options += ' --max -E 1000 --incE 1000 '
            self.binary = 'hmmscan'  # protein version (hmmscan rather than nhmmscan) seems to be more sensitive
        else:
            assert self.sensitivity == ''

        std_out_treatment = ''
        if self.binary == 'nhmmscan':
            options += ' --tblout >(sed s/^/-TBLOUT-/)' #' --tblout >(cat)'
            std_out_treatment = '' #'>/dev/null'  # for nhmmscan, we only need the tblout output, so redirect throw out the regular output
        elif self.binary == 'hmmscan':  # whereas for some **!#&$! reason hmmscan doesn't report the alignment start and end positions in its tfjfjfjffj
            options += ' --domtblout >(sed s/^/-TBLOUT-/)'
            std_out_treatment = ''
        else:
            assert False  # er, case not covered, so yer prolly screwed anyway
        # here the '-' is the seq file:
        return query_seq_command + ' | /bin/bash -c \"' + hmmerdir + '/binaries/' + self.binary + ' ' + options + ' ' + hmm_dbase_file + ' - ' + std_out_treatment + '\"'

    def run(self):
#        print 'running hmmer with %s' % self.seq
        self.output = check_output(self.command, shell=True)
        self.parse_output()
#        if self.sensitivity == 'max' and len(self.matches) == 0:
#            print '\nno matches with max sensitivity\n'
#            check_call(self.command.replace('>/dev/null','| grep IGH'), shell=True)
#            sys.exit()

    def parse_output(self):
        assert len(self.matches) == 0  # um, no particular reason to suspect otherwise a.t.m.
        outlines = self.output.splitlines()
        outlines = dropwhile(lambda line: line.find('-TBLOUT-') < 0, outlines)
        header_line = sanitize_header_line(outlines.next().replace('-TBLOUT-',''))
            
        if self.binary == 'hmmscan':
#            outlines.next()
            header_line = sanitize_header_line(outlines.next().replace('-TBLOUT-',''))  # NOTE this is the overall *sequence* evalue I'm pulling out here. I don't care a.t.m., but keep it in mind
            header_line = sanitize_hmmscan_header_line(header_line)
        # make a dict so we don't have to remember the order
        indices = {}
        for icol in range(len(header_line.split())):
            column = header_line.split()[icol]
            indices[column] = icol
        for outline in outlines:
            outline = outline.replace('-TBLOUT-','')
            if outline[0] == '#':
                continue
            outline = outline.split()
            entry = {}
            entry['target_name'] = outline[indices['target_name']]
            for column in ['hmm_from', 'hmm_to', 'ali_from', 'ali_to']:
                entry[column] = int(outline[indices[column]])
            entry['evalue'] = float(outline[indices['evalue']])
            if entry['ali_from'] > entry['ali_to']:
                print '    skipping reverse match %d %d' % (entry['ali_from'], entry['ali_to'])
                continue
            self.matches.append(entry)
            if len(self.matches) > 1:  # only look at the first two matches a.t.m.
                break

""" Execute search with proper hmmer binary, and store output in a dict. """

from subprocess import check_output

#----------------------------------------------------------------------------------------
def sanitize_header_line(line):
    """ Really, dude? You couldn't make the damn headers one word each? """
    line = line.replace('#', '')
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
class Hmmerer(object):
    """ ya what it says up there yo """

    def __init__(self, hmmerdir, region, seq, sensitivity=''):
        self.region = region
        self.seq = seq
        self.sensitivity = sensitivity
        self.command = self.build_command(hmmerdir)
        self.output = ''  # raw output from hmmer
        self.matches = []

    def build_command(self, hmmerdir):
        """ Build hmmer command. """
        hmm_dbase_file = 'data/hmms/' + self.region + '/all.hmm'
        query_seq_command = 'echo \"> test seq file\n' + self.seq + '\n\"'  # command to generate query sequence 'file': pipe it to hmmer to avoid writing to disk
        binary = '/bin/bash -c \"' + hmmerdir + '/binaries/nhmmscan'
        options = '--tblout >(cat)'
        if self.sensitivity == 'max':
            options += ' --max'
            # binary = binary.replace('nhmmscan', 'hmmscan')  # docs are a little unclear on the real differences -- but in practice hmmscan gets more matches
        else:
            assert self.sensitivity == ''

        return query_seq_command + ' | ' + binary + ' ' + options + ' ' + hmm_dbase_file + ' - >/dev/null\"'

    def run(self):
        self.output = check_output(self.command, shell=True)
        self.parse_output()

    def parse_output(self):
        assert len(self.matches) == 0  # um, no particular reason to suspect otherwise a.t.m.
        outlines = self.output.splitlines()
        header_line = sanitize_header_line(outlines[0])
        # make a dict so we don't have to remember the order
        indices = {}
        for icol in range(len(header_line.split())):
            column = header_line.split()[icol]
            indices[column] = icol
        for outline in outlines:
            if outline[0] == '#':
                continue
            outline = outline.split()
            entry = {}
            entry['target_name'] = outline[indices['target_name']]
            for column in ['hmm_from', 'hmm_to', 'ali_from', 'ali_to']:
                entry[column] = int(outline[indices[column]])
            entry['evalue'] = float(outline[indices['evalue']])
            self.matches.append(entry)
            if len(self.matches) > 1:  # only look at the first two matches a.t.m.
                break
    
    #    if len(matches) == 0:  # no matches
    #        print ' (none)',
    #    else:
    #        print ' (' + matches[0]['target_name']+ ')',
            
            

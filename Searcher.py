""" Run hmmer via the Hmmerer class, and make all the decisions about search order, strategies, heuristics, & co. """

import sys
from itertools import takewhile, dropwhile, tee

from Hmmerer import Hmmerer
import utils

class Searcher(object):
    """ same ol shit """
    def __init__(self, initial_query_seq, debug=False, n_matches_max=999):
        self.n_matches_max = n_matches_max
        self.debug = debug
        self.hmmerdir = '/home/dralph/Dropbox/work/hmmer/hmmer-3.1b1-linux-intel-x86_64'
        self.seq = initial_query_seq
        self.query_seqs = {}  # the seq that was actually used for each query (i.e. without excised portions)
        self.current_seq = self.seq  # we will in general excise parts of the original sequence, and store the result in current_seq
        self.n_tries = {}  # how many tries for this region?
        self.matches = {}
        self.best_matches = {}
        self.excisions = []
        for region in utils.regions:
            self.n_tries[region] = 0
            self.query_seqs[region] = ''
            self.matches[region] = []
            self.best_matches[region] = {}
        
    def all_matched(self):
        """ Do we have at least one match for each region? """
        return ((len(self.best_matches['v']) > 0) and
                (len(self.best_matches['d']) > 0) and
                (len(self.best_matches['j']) > 0))
        
    def search(self):
        found_str = ''
#        print '.',
        if self.debug:
            print '\nsearching'
        for region in ['v', 'j', 'd']:  # do d last
            self.get_matches(region)
            if len(self.matches[region]) != 0 and len(self.best_matches[region]) != 0:
                found_str += region
                self.excise_match(region)   # if there's a good match cut it out. TODO change 'good' criterion
                if len(self.current_seq) == 0:
                    return
        # now try again
        for region in ['v', 'j', 'd']:  # do d last, still
            if not self.is_matched(region):
                self.get_matches(region)
                if len(self.matches[region]) != 0 and len(self.best_matches[region]) != 0:
                    found_str += region
                    self.excise_match(region)   # if there's a good match cut it out. TODO change 'good' criterion
                    if len(self.current_seq) == 0:  # excised the whole sequence
                        return
                else:
                    pass  #print '  WARNING no matches found',
        return found_str

    def get_matches(self, region, debug=False):
        """ Look for matches for <region>. """
        sensitivity = ''
        if self.n_tries[region] == 1:
            sensitivity = 'more'
        elif self.n_tries[region] > 1:
            sensitivity = 'max'
        hmmerer = Hmmerer(self.hmmerdir, region, self.current_seq, sensitivity=sensitivity, debug=self.debug, n_matches_max=self.n_matches_max)
        hmmerer.run()
#        print hmmerer.output
        self.n_tries[region] += 1
        self.matches[region] = hmmerer.matches
        self.find_best_match(hmmerer, region)

    def is_matched(self, region):
        """ Have we found matches for this region? """
        return len(self.matches[region]) != 0

    def get_best_match_name(self, region):
        if len(self.matches[region]) == 0 or len(self.best_matches[region]) == 0:
            return 'none'
        else:
            return self.best_matches[region]['target_name']

    def find_best_match(self, hmmerer, region):
        if self.debug:
            print '  looking for best match of %d' % len(self.matches[region])
        best_evalue = 999.0
        for match in self.matches[region]:
            if match['ali_from'] > match['ali_to']:  # NOTE these are index *one* counting (!!!)
                print '  REVERSE %d --> %d (skipping)' % (match['ali_from'],match['ali_to'])
                continue
            self.get_matching_section(hmmerer, match)  # fill output section in this match
            if '.' in match['hmm_seq']:  # remove matches with deletions from consideration
                if self.debug:
                    print '  eliminating match with deletion %s' % match['hmm_seq'].upper()
                continue
            if '-' in match['test_seq']:  # remove matches with deletions from consideration
                if self.debug:
                    print '  eliminating match with deletion %s' % match['test_seq']
                continue
            if match['evalue'] < best_evalue:
                best_evalue = match['evalue']
                self.best_matches[region] = match

    def excise_match(self, region):
        """ Remove best_match from current_seq.

        Once we're sure enough we have the correct germ line gene
        for one region {vdj}, cut it out to help with identification
        of the others.
        """
    
        # convert from 1-start to 0-start indexing (!!!)
        excise_from = self.best_matches[region]['ali_from'] - 1
        excise_to = self.best_matches[region]['ali_to'] - 1

        if region == 'v' and excise_from != 0:
            excise_from = 0
            assert len(self.current_seq[:excise_from]) == 0
            if self.debug:
                print '   expanding v excision to zero'
        if region == 'j' and excise_to != len(self.current_seq) - 1:
            excise_to = len(self.current_seq) - 1
            assert len(self.current_seq[excise_to + 1 :]) == 0
            if self.debug:
                print '   expanding j excision to %d' % (len(self.current_seq) - 1)

        self.excisions.append({'region': region, 'from': excise_from, 'to': excise_to})
        if self.debug:
            print '    excising from %d to %d: %s --> %s' % (excise_from, excise_to, self.current_seq, self.current_seq[:excise_from] + self.current_seq[excise_to + 1 :])

        self.query_seqs[region] = self.current_seq
        self.current_seq = self.current_seq[:excise_from] + self.current_seq[excise_to + 1 :]

    def build_inferred_seq(self, seq, all_germlines, outline):
        assert self.excisions[0]['region'] == 'v'  # makes it easier a.t.m.
        assert self.excisions[1]['region'] == 'j'
        assert self.excisions[2]['region'] == 'd'
        germlines, hmms, ihmms = {}, {}, {}
        for region in utils.regions:
            germlines[region] = all_germlines[region][utils.unsanitize_name(self.best_matches[region]['target_name'])]
            hmms[region] = self.best_matches[region]['hmm_seq']
            ihmms[region] = germlines[region].find(hmms[region].upper())  # position at which the consensus (hmm) starts in the germline sequence
            try:
                assert ihmms[region] >= 0
            except:
                print germlines[region]
                print hmms[region].upper()
                print ihmms[region]
                assert False
            print '  hmm for %s runs from %d to %d (inclusive)' % (region, ihmms[region], ihmms[region] + len(hmms[region]) - 1)

        outline['v_5p_del'] = ihmms['v'] # TODO kinda otter be zero
        outline['v_3p_del'] = len(germlines['v']) - ihmms['v'] - len(hmms['v']) # len(germlines['v']) - len(hmms['v']) - germlines['v'].find(hmms['v'].upper())
        outline['d_5p_del'] = ihmms['d'] # germlines['d'].find(hmms['d'].upper())
        outline['d_3p_del'] = len(germlines['d']) - ihmms['d'] - len(hmms['d']) # len(germlines['d']) - len(hmms['d']) - germlines['d'].find(hmms['d'].upper())
        outline['j_5p_del'] = ihmms['j'] # germlines['j'].find(hmms['j'].upper())
        outline['j_3p_del'] = len(germlines['j']) - ihmms['j'] - len(hmms['j'])  # TODO kinda otter be zero

        for ex in self.excisions:
            match = self.best_matches[ex['region']]
            print '  excised match %s: %d --> %d' % (ex['region'], ex['from'], ex['to'])
            print '        test %s' % match['test_seq']
            hmm_start = ihmms[ex['region']]
            hmm_end = ihmms[ex['region']] + len(hmms[ex['region']]) - 1
            print '         hmm %s' % (hmm_start * '.' + match['hmm_seq'].upper() + (len(germlines[ex['region']]) - ihmms[ex['region']] - len(hmms[ex['region']])) * '.')  # NOTE ali_from includes the d part!
            print '    germline %s' % all_germlines[ex['region']][utils.unsanitize_name(match['target_name'])]

        #----------------------------------------------------------------------------------------
        # NOTE these are inclusive
        seq_match_start, seq_match_end = {}, {}
        seq_match_start['v'] = self.best_matches['v']['ali_from'] - 1
        seq_match_end['v'] = seq_match_start['v'] + len(hmms['v']) - 1
        seq_match_start['d'] = self.excisions[0]['to'] - self.excisions[0]['from'] + self.best_matches['d']['ali_from']
        seq_match_end['d'] = seq_match_start['d'] + len(hmms['d']) - 1
        seq_match_start['j'] = self.excisions[0]['to'] - self.excisions[0]['from'] + self.best_matches['j']['ali_from']
        seq_match_end['j'] = seq_match_start['j'] + len(hmms['j']) - 1
        outline['vd_insertion'] = seq[seq_match_end['v']+1 : seq_match_start['d']]
        outline['dj_insertion'] = seq[seq_match_end['d']+1 : seq_match_start['j']]

        actual_seq_length = len(seq)
        inferred_seq_length = outline['v_5p_del'] + len(hmms['v']) + len(outline['vd_insertion']) + len(hmms['d']) + len(outline['dj_insertion']) + len(hmms['j']) + outline['j_3p_del']
        print '    actual %d  inferred %d' % (actual_seq_length,inferred_seq_length)
        if actual_seq_length != inferred_seq_length:
            
            outline['ack'] = True

    def get_matching_section(self, hmmerer, match):
        if len(match) == 0:
            return
        if self.debug:
            print '  getting detailed info for best match %s' % match['target_name']
        output = dropwhile(lambda line: line.find('>> ' + match['target_name']) < 0, hmmerer.output.splitlines())  # drop until you get to the block on this gene
        try:
            output.next()
        except:
            print hmmerer.output
            assert False
        output = takewhile(lambda line: line.find('>> IGH') < 0, output)  # drop the next block
        if self.debug:
            output, debug_output = tee(output)
            for line in debug_output:
                print '     -->    ',line
        output = dropwhile(lambda line: line.find(match['target_name']) < 0, output)  # drop the header crap
        output = takewhile(lambda line: line.find('== domain') < 0, output)  # drop domains after the first one. TODO don't do that
        output = takewhile(lambda line: line.find('Internal pipeline') < 0, output)  # kill pipeline statistics
        hmm_seq, matching_position_line, test_seq, probabilities = '', '', '', ''
        for line in output:
            # TODO I think this will fail if there's indels in the interior
            if len(line) == 0:
                continue
            # first get the germline line
            line = line.split()
            assert line[0] == match['target_name']
            start_position = line[1]  # start position of this match
            hmm_seq += line[2]
            end_position = line[3]
            # then get the line with positions at which they match
            matching_position_line +=  output.next()  # screwed up with extra spaces at the start, but I'm not using it a.t.m.
            # then the line for the test sequence
            line = output.next().split()
            assert line[0] == 'test'
#            assert line[1] == start_position  # um, maybe not?
            test_seq += line[2]
#            assert line[3] == end_position  # um, maybe not?
            # and finally, get the probability line
            probabilities += output.next().split()[0]
        match['hmm_seq'] = hmm_seq
#        match[region + '_5p_del'] = match['hmm_from'] - 1
#        # TODO is hmm the same length as the germline gene?
#        match[region + '_3p_del'] = 0 #len(hmm_seq) - match['hmm_to'] - 1
        match['test_seq'] = test_seq

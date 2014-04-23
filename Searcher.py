""" Run hmmer via the Hmmerer class, and make all the decisions about search order, strategies, heuristics, & co. """

import partutils
from Hmmerer import Hmmerer

class Searcher(object):
    """ same ol shit """
    def __init__(self, query_line):
        self.hmmerdir = '/home/dralph/Dropbox/work/hmmer/hmmer-3.1b1-linux-intel-x86_64'
        self.seq = query_line['seq']
        self.current_seq = self.seq  # we will in general excise parts of the original sequence, and store the result in current_seq
        self.n_tries = {}  # how many tries for this region?
        self.matches = {}
        self.best_matches = {}
        for region in partutils.regions:
            self.n_tries[region] = 0
            self.matches[region] = []
            self.best_matches[region] = {}
        
    def search(self):
        for region in ['v', 'j', 'd']:  # do d last
            self.get_matches(region)
            if len(self.matches[region]) != 0:
                self.excise_match(region)   # if there's a good match cut it out. TODO change 'good' criterion
        # now try again!
        for region in ['v', 'j', 'd']:  # do d last, still
            if not self.is_matched(region):
                self.get_matches(region)
                if len(self.matches[region]) != 0:
                    self.excise_match(region)   # if there's a good match cut it out. TODO change 'good' criterion
            if not self.is_matched(region):
                print '\n  dammit %s still not matched: %s' % (region, self.current_seq),
#                get_matches(hmmdir, region, current_seq, is_matched, n_tries, inline[region + '_gene'], is_correct, debug=True)

    def get_matches(self, region, debug=False):
        """ Look for matches for <region>. """
        sensitivity = ''
        if self.n_tries[region] > 0:
            sensitivity = 'max'
        hmmerer = Hmmerer(self.hmmerdir, region, self.current_seq, sensitivity=sensitivity)
        hmmerer.run()
        self.n_tries[region] += 1
        self.matches[region] = hmmerer.matches
        self.find_best_match(region)

    def is_matched(self, region):
        """ Have we found matches for this region? """
        return len(self.matches[region]) != 0

    def get_best_match_name(self, region):
        if len(self.matches[region]) == 0:
            return 'none'
        else:
            return self.best_matches[region]['target_name']

    def find_best_match(self, region):
        if len(self.matches[region]) == 1:
            self.best_matches[region] = self.matches[region][0]
        else:
            best_evalue = 999.0
            for match in self.matches[region]:
                if match['evalue'] < best_evalue:
                    best_evalue = match['evalue']
                    self.best_matches[region] = match

        if len(self.matches[region]) > 0:
               len(self.best_matches[region]) != 0

    def excise_match(self, region):
        """ Remove best_match from current_seq.

        Once we're sure enough we have the correct germ line gene
        for one region {vdj}, cut it out to help with identification
        of the others.
        """
    
        # convert from 1-start to 0-start indexing
        excise_from = self.best_matches[region]['ali_from'] - 1
        excise_to = self.best_matches[region]['ali_to'] - 1

        self.current_seq = self.current_seq[:excise_from] + self.current_seq[excise_to + 1 :]

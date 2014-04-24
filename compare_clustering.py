#!/usr/bin/env python
"""
Compare two clusterings.
Each input should be a CSV file with at least columns containing a sequence
name and (integer) group identifier.
"""
import argparse
import collections
import csv
import logging
import json
import sys
import os
import bz2

import numpy as np
from sklearn import metrics


NAME_DEFAULT = 'name'
GROUP_DEFAULT = 'group'
_EXT_MAP = {'.bz2': bz2.BZ2File}

_log = logging.getLogger('compare-clustering')


def handle_factory(mode='r'):
    def opener(s):
        if s == '-' and mode == 'r':
            return sys.stdin
        elif s == '-' and mode == 'w':
            return sys.stdout
        e = os.path.splitext(s)[1]
        return _EXT_MAP.get(e, open)(s, mode=mode)
    return opener


def load_clustering(fp, name_column=NAME_DEFAULT, group_column=GROUP_DEFAULT,
                    has_name=True, **kwargs):
    """Load a clustering from a file"""
    reader = csv.reader(fp, **kwargs)
    col_names = next(reader, None)
    result = {}
    group_index = col_names.index(group_column)
    if not has_name:
        for i, row in enumerate(reader):
            result[str(i)] = long(row[group_index])
    else:
        name_index = col_names.index(name_column)

        for row in reader:
            result[row[name_index]] = long(row[group_index])
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument('true_csv', type=handle_factory('r'),
                   help="""CSV file containing true clustering.""")
    p.add_argument('inferred_csv', type=handle_factory('r'),
                   help="""CSV file containing inferred clusters""")

    true_group = p.add_argument_group('True CSV file')
    true_group.add_argument('--true-name-column', default=NAME_DEFAULT,
                            help="""Name of sequence""")
    true_group.add_argument('--true-group-column', default=GROUP_DEFAULT,
                            help="""Name of sequence""")
    true_group.add_argument('--no-name', default=True, action='store_false',
                            dest='has_name', help="""Sequences do not have
                            names - label numerically starting at 0.""")
    inferred_group = p.add_argument_group('Inferred CSV file')
    inferred_group.add_argument('--inferred-name-column', default=NAME_DEFAULT,
                                help="""Name of sequence""")
    inferred_group.add_argument('--inferred-group-column',
                                default=GROUP_DEFAULT, help="""Name of
                                sequence""")

    out_grp = p.add_argument_group('Output')
    out_grp.add_argument('-o', '--outfile', type=handle_factory('w'),
                         default=sys.stdout)
    out_grp.add_argument('-f', '--format', default='json', help="""Output
                         format""", choices=('csv', 'json'))

    a = p.parse_args()

    with a.true_csv as fp:
        true_clustering = load_clustering(fp, name_column=a.true_name_column,
                                          group_column=a.true_group_column,
                                          has_name=a.has_name)
    with a.inferred_csv as fp:
        inferred_clustering = load_clustering(
            fp,
            name_column=a.inferred_name_column,
            group_column=a.inferred_group_column)

    singletons = frozenset(true_clustering) ^ frozenset(inferred_clustering)
    if singletons:
        raise ValueError(
            'Non-empty symmetric difference: {0}'.format(singletons))

    names = sorted(true_clustering)
    true_arr = np.array([true_clustering[name] for name in names])
    inferred_arr = np.array([inferred_clustering[name] for name in names])
    _log.info('Read clusters. Calculating MI.')

    result = [('mi', metrics.mutual_info_score(true_arr, inferred_arr)),
              ('normalized_mi',
               metrics.normalized_mutual_info_score(true_arr, inferred_arr)),
              ('adjusted_mi',
               metrics.adjusted_mutual_info_score(true_arr, inferred_arr)),
              ('f1_score',
               metrics.f1_score(true_arr, inferred_arr)),
              ('completeness_score',
               metrics.completeness_score(true_arr, inferred_arr)),
              ('homogeneity_score',
               metrics.homogeneity_score(true_arr, inferred_arr))]

    with a.outfile as ofp:
        if a.format == 'csv':
            writer = csv.writer(ofp, lineterminator='\n')
            writer.writerow(['metric', 'value'])
            writer.writerows(result)
        elif a.format == 'json':
            odict = collections.OrderedDict
            rec = odict()
            rec['metrics'] = odict(result)
            rec['true_labels'] = list(true_arr)
            rec['inferred_labels'] = list(inferred_arr)
            json.dump(rec, ofp, indent=2)
            ofp.write('\n')
        else:
            raise ValueError('Unknown format: {0}'.format(a.format))


if __name__ == '__main__':
    main()

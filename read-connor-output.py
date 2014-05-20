#!/usr/bin/env python
import sys
import os
import bz2
import csv
from opener import opener
import utils

data_type = sys.argv[1]
#infname = 'head-' + data_type + '.csv'  # bzgrep -m100 . /shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/C/01-C-N_merged.tsv.bz2 | sed 's/[ \t][ \t]*/,/g'|cut -f2 -d, |sed 's/nucleotide/seq/'> head-data.csv

infname = '/shared/silo_researcher/Matsen_F/MatsenGrp/working/cmccoy/20140414-bcell-validation/build/out/annotations.csv.bz2'

n_ties = {}
for region in utils.regions:
    n_ties[region] = []

with opener('r')(infname) as infile:
    reader = csv.DictReader(infile)
    for line in reader:
        for region in utils.regions:
            tie_column = region.upper() + '_ties'
            number = 0
            if len(line[tie_column]) > 0:
                number = len(line[tie_column].split('|'))
            n_ties[region].append(number)

for region in utils.regions:
    var = 'ties'
    outdir = 'data/' + data_type + '/' + var
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    with opener('w')(outdir + '/' + region+ '.txt') as outfile:
        for value in n_ties[region]:
            outfile.write('%7e\n' % value)

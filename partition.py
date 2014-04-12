#!/usr/bin/env python
import csv
from subprocess import call
from glob import glob
import os
import sys

from opener import opener

regions = ('v','d','j')

def build_hmms():
    hmmdir = '/home/dralph/work/recombinator/data/msa'
    for region in regions:
        for sto_fname in glob(hmmdir + '/' + region + '/*.sto'):
#            with opener('r')(sto_fname) as sto_file:
            outdir = 'data/hmms/' + region
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            call(["../binaries/hmmbuild", '--dna', outdir + '/' + os.path.basename(sto_fname).replace('.sto','.hmm'), sto_fname]) #, stdout=foo)

#seq_data = []
#with open('/home/dralph/work/recombinator/out.csv') as infile:
#    reader = csv.reader(infile)
#    for row in reader:
#        seq_data.append({'hash':row[0], 'v_gene':row[2], 'd_gene':row[3], 'j_gene':row[4], 'seq':row[13]})
#
#for row in seq_data:
#    if row['seq'] == 'seq':
#        continue
#    print row['hash'],row['seq']
#    # NOTE I shouldn't have to make a stupid temp file, but I can't get it to accept stdin
#    with open('/tmp/partition-tmp.cfg', 'w') as cfg_file:
#        cfg_file.write('> ' + row['hash'] + '\n')
#        cfg_file.write(row['seq'] + '\n')
#    foo = ''
#    call(["../binaries/hmmscan", 'all-vdjs.hmm', '/tmp/partition-tmp.cfg'], stdout=foo)
#    assert False

build_hmms()

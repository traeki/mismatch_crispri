#!/usr/bin/env python

# Author: John Hawkins (jsh) [really@gmail.com]

import argparse
import collections
import gzip
import logging
import random
import sys

from Bio import SeqIO
from Bio import Seq


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')


def parse_args():
  logging.info('Parsing command line.')
  parser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--guide_set', type=str, required=True,
                      help='Location of intended guide list.')
  parser.add_argument('--input_fastq', type=str, required=True,
                      help='Location of read file in FASTQ format.')
  parser.add_argument('--reverse', action='store_true',
                      help='If set, guide is oriented opposite to read direction.')
  args = parser.parse_args()
  return args


def main():
  args = parse_args()
  if args.input_fastq.endswith('.gz'):
    handle = gzip.open(args.input_fastq, 'rt')
  else:
    handle = open(args.input_fastq, 'r')
  hitlist = set([x.strip() for x in open(args.guide_set, 'r')])
  outfile = open(args.input_fastq + '.counts', 'w')
  weirdfile = open(args.input_fastq + '.weird', 'w')
  skipfile = open(args.input_fastq + '.skipped', 'w')
  skipped = list()
  counts = collections.defaultdict(int)
  reads = 0
  for x in hitlist:
    counts[x] = 0
  for i, record in enumerate(SeqIO.parse(handle, 'fastq-sanger')):
    reads += 1
    if random.random() < 0.00001:
      logging.info('considering record {0}: {1}'.format(i, record))
    s = record.seq
    if not args.reverse:
      endpos = s.find('GTTTTAGAG')
      startpos = 0
    else:
      header = 'TCTAAAAC'
      startpos = s.find(header)
      if startpos >= 0:
        startpos += len(header)
        endpos = startpos + 20
    if startpos >= 0 and endpos >= 0:
      s = str(s[startpos:endpos])
      if args.reverse:
        s = str(Seq.Seq(s).reverse_complement())
      counts[s] += 1
    else:
      skipped.append(record)
    if len(skipped) > 10000:
      logging.info('DUMPING SKIPPED RECORDS')
      SeqIO.write(skipped, skipfile, 'fastq-sanger')
      skipped = list()
  logging.info('DUMPING FINAL SKIPPED RECORDS')
  SeqIO.write(skipped, skipfile, 'fastq-sanger')
  logging.info('Sorting records')
  for k,v in sorted(iter(counts.items()), key=lambda k_v: k_v[1], reverse=True):
    if k in hitlist:
      outfile.write('\t'.join([k, str(v)]) + '\n')
    else:
      weirdfile.write('\t'.join([k, str(v)]) + '\n')
  handle.close()
  hits = len(hitlist)
  if hits == 0:
    ratio == 'n/a'
  else:
    ratio = reads/hits
  logging.info('mean(reads/oligo) = {0}/{1} = {2}'.format(reads, hits, ratio))

##############################################
if __name__ == "__main__":
  sys.exit(main())

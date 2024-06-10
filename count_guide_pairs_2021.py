#!/usr/bin/env python

# Author: John Hawkins (jsh) [really@gmail.com]

import argparse
import collections
import gzip
import itertools
import logging
import os.path
import random
import string
import sys

from Bio import SeqIO
from Bio import Seq


# logging.basicConfig(level=logging.DEBUG,
#                     format='%(asctime)s %(levelname)s %(message)s')
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

class Error(Exception):
  pass

class SampleError(Error):
  pass


def parse_locus_map(map_file):
  locus_map = dict()
  for line in open(map_file, 'r'):
    guide, bc, locus = [x.strip() for x in line.split('\t')]
    locus_map[bc[1:]] = locus
    locus_map[guide] = locus
  return locus_map


def parse_args():
  """Read in the arguments for the sgrna library construction code."""
  logging.info('Parsing command line.')
  parser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('--locus_map', type=str, required=True,
                      help='Location of guide/barcode mapping.')
  parser.add_argument('--front_fastq', type=str, required=True,
                      help='Location of front read file in FASTQ format.')
  parser.add_argument('--rear_fastq', type=str, required=True,
                      help='Location of rear read file in FASTQ format.')
  args = parser.parse_args()
  # if args.tsv_file_name is None:
  #   base = os.path.splitext(args.input_fasta_genome_name)[0]
  #   args.tsv_file_name =  base + '.targets.all.tsv'
  return args


def main():
  args = parse_args()
  front_handle = gzip.open(args.front_fastq, 'rt')
  front_records = SeqIO.parse(front_handle, 'fastq-sanger')
  rear_handle = gzip.open(args.rear_fastq, 'rt')
  rear_records = SeqIO.parse(rear_handle, 'fastq-sanger')
  outfile = open(args.front_fastq + '.counts', 'w')
  pairfile = open(args.front_fastq + '.pairs', 'w')
  frontfile = open(args.front_fastq + '.front', 'w')
  rearfile = open(args.front_fastq + '.rear', 'w')
  weirdfile = open(args.front_fastq + '.weird', 'w')
  skip_front = open(args.front_fastq + '.skipped', 'w')
  skip_rear = open(args.rear_fastq + '.skipped', 'w')
  skipped_front = list()
  skipped_rear = list()
  counts = collections.defaultdict(int)
  sample_fraction = 1./1000
  record_i = 0
  for front, rear in zip(front_records, rear_records):
    if random.random() < sample_fraction:
      id = front.id
      logging.info('processing record {record_i}: {id}'.format(**vars()))
    record_i += 1
    assert front.id == rear.id
    # find the first guide (from read2)
    f = front.seq
    f_startpos = 1 # NOTE(jsh): first base is always N, locus_map modified to match
    f_endpos = f.find('ATAGGGAACT')
    # find the second guide (from read1)
    r = rear.seq
    r_header = 'GCTCTTAAAC'
    r_startpos = r.find(r_header)
    if r_startpos >= 0:
      r_startpos += len(r_header)
    # r_endpos = r.find('ACATAGATTA') <-- for old pre-BMK XY library
    r_endpos = r.find('ACATTAAGTA')
    if f_endpos < 0 or f_startpos < 0 or r_endpos < 0 or r_startpos < 0:
      skipped_front.append(front)
      skipped_rear.append(rear)
    else:
      f = str(f[f_startpos:f_endpos])
      r = str(r[r_startpos:r_endpos])
      r = str(Seq.Seq(r).reverse_complement())
      counts[(f, r)] += 1
    if len(skipped_front) > 10000:
      logging.info('DUMPING SKIPPED RECORDS')
      SeqIO.write(skipped_front, skip_front, 'fastq-sanger')
      SeqIO.write(skipped_rear, skip_rear, 'fastq-sanger')
      skipped_front = list()
      skipped_rear = list()
  logging.info('DUMPING FINAL SKIPPED RECORDS')
  SeqIO.write(skipped_front, skip_front, 'fastq-sanger')
  SeqIO.write(skipped_rear, skip_rear, 'fastq-sanger')
  skipped_front = list()
  skipped_rear = list()
  # set up locus_map / expected
  expected = dict()
  locus_map = parse_locus_map(args.locus_map)
  for x in list(locus_map.values()):
    for y in list(locus_map.values()):
      expected[(x,y)] = 0
  front_stats = dict([(x, 0) for x in list(locus_map.values())])
  rear_stats = dict([(x, 0) for x in list(locus_map.values())])
  pair_stats = collections.defaultdict(int)
  for k,v in sorted(iter(list(counts.items())), key=lambda k_v: k_v[1], reverse=True):
    f, r = k
    if r not in locus_map or f not in locus_map:
      weirdfile.write('\t'.join([f, r, str(v)]) + '\n')
      continue
    f = locus_map[f]
    r = locus_map[r]
    front_stats[f] += v
    rear_stats[r] += v
    pair = f <= r and (f,r) or (r,f)
    pair_stats[pair] += v
    expected[(f, r)] = v
  for pair, count in list(expected.items()):
    a, b = pair
    outfile.write('\t'.join((a, b, str(count))) + '\n')
  for pair, count in list(pair_stats.items()):
    a, b = pair
    pairfile.write('\t'.join((a, b, str(count))) + '\n')
  for front, count in list(front_stats.items()):
    frontfile.write('\t'.join((front, str(count))) + '\n')
  for rear, count in list(rear_stats.items()):
    rearfile.write('\t'.join((rear, str(count))) + '\n')
  # output a tsv
  SeqIO.write(skipped_front, skip_front, 'fastq-sanger')
  SeqIO.write(skipped_rear, skip_rear, 'fastq-sanger')

##############################################
if __name__ == "__main__":
  sys.exit(main())

#!/usr/bin/env python
# Author: John Hawkins (jsh) [really@gmail.com]

import argparse
import logging
import pathlib
import shutil
import sys

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as st
import seaborn as sns

import eval_lib as el

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s') 

_PACKAGEDIR = pathlib.Path(__file__).parent
TESTDIR = _PACKAGEDIR / 'testdata'
_CODEFILE = pathlib.Path(__file__)


def parse_args():
  """Read in the arguments for the sgrna library construction code."""
  logging.info('Parsing command line.')
  parser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument(
      '--gammas', type=str,
      help='file: unique entries for gamma by variant -- x axis',
      required=True)
  parser.add_argument(
      '--name', type=str,
      help='str: label name for x axis',
      required=True)
  parser.add_argument(
      '--maxguides', type=int,
      help='int: largest number of guides to count in bin',
      required=True)
  parser.add_argument(
      '--pngfile', type=str,
      help='file: output name for figure',
      required=True)
  args = parser.parse_args()
  return args


def main():
  args = parse_args()

  template = 'Reading {args.name} from {args.gammas}...'
  logging.info(template.format(**locals()))
  data = pd.read_csv(args.gammas, sep='\t')
  bintally = el.bin_counts_by_gene(data)


  fig, ax = plt.subplots(1, 1, figsize=(6,3), sharex=True)
  gbin = bintally['mid'].fillna(0)
  nbins = min(args.maxguides, 20)
  bins = np.linspace(0,args.maxguides,nbins-1)
  sns.distplot(gbin, ax=ax, bins=bins, kde=False)
  binmed = int(gbin.median())
  binlow = int(gbin.quantile(0.05))
  ax.text(0.05, 0.80,
          'median: {binmed}\n95% >= {binlow}'.format(**locals()),
          transform=ax.transAxes)
  plt.xlim(0,args.maxguides)
  plt.xlabel('# of guides w/ 10% to 90% wildtype growth rate')
  plt.ylabel('# of genes')
  plt.tight_layout()
  plt.savefig(args.pngfile, dpi=600)
  plt.close('all')

##############################################
if __name__ == "__main__":
  sys.exit(main())

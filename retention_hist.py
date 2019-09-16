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
import gamma_lib as gl

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
      '--startcounts', type=str,
      help='file: unique entries for start count by variant',
      required=True)
  parser.add_argument(
      '--annofile', type=str,
      help='file: gamma file for gene annotations',
      required=True)
  parser.add_argument(
      '--pngfile', type=str,
      help='file: output name for figure',
      required=True)
  args = parser.parse_args()
  return args


def main():
  args = parse_args()

  template = 'Reading counts from {args.startcounts}...'
  logging.info(template.format(**locals()))
  data = pd.read_csv(args.startcounts, sep='\t', header=None)
  data.columns = ['variant', 'count']
  anno = pd.read_csv(args.annofile, sep='\t')
  data = data.merge(anno[['variant', 'gene']], on='variant', how='left')
  threshold = gl.MIN_START_READS
  data['overline'] = (data['count'] > threshold)
  genevals = data.groupby('gene').overline.sum()
  pergene_file = pathlib.Path(args.pngfile)
  overall_file = pergene_file.with_suffix('.overall.png')
  data['logcount'] = np.log2(data['count'].clip(1))

  import IPython; IPython.embed()
  sys.exit(1)

  fig, ax = plt.subplots(1, 1, figsize=(6,3), sharex=True)
  bins = np.linspace(0,100,101)
  sns.distplot(genevals, ax=ax, bins=bins, kde=False)
  binmed = int(genevals.median())
  binlow = int(genevals.quantile(0.05))
  ax.text(0.45, 0.05,
          'median: {binmed}\n95% >= {binlow}'.format(**locals()),
          transform=ax.transAxes)
  template = 'number of guides with start count > {threshold}'
  plt.xlabel(template.format(**locals()))
  plt.ylabel('number of genes in bin')
  plt.tight_layout()
  plt.savefig(pergene_file, dpi=600)
  plt.close('all')

  fig, ax = plt.subplots(1, 1, figsize=(6,3), sharex=True)
  sns.distplot(data.logcount, ax=ax, kde=False)
  plt.plot([np.log2(threshold), np.log2(threshold)], [0, 2000], 'r-')
  template = 'raw count of guide at t0 (log2)'
  plt.xlabel(template.format(**locals()))
  plt.ylabel('number of guides in bin')
  plt.tight_layout()
  plt.savefig(overall_file, dpi=600)
  plt.close('all')


##############################################
if __name__ == "__main__":
  sys.exit(main())

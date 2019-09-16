#!/usr/bin/env python
# Author: John Hawkins (jsh) [really@gmail.com]

import argparse
import logging
import pandas as pd
import pathlib
import shutil
import sys

from matplotlib import pyplot as plt
import seaborn as sns
import scipy.stats as st

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s') 

_PACKAGEDIR = pathlib.Path(__file__).parent
TESTDIR = _PACKAGEDIR / 'testdata'
_CODEFILE = pathlib.Path(__file__)


def parse_args():
  logging.info('Parsing command line.')
  parser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument(
      '--gammas', type=str,
      help='file: unaveraged gammas by replicate',
      required=True)
  parser.add_argument(
      '--pngfile', type=str,
      help='file: name of file for plot',
      required=True)
  args = parser.parse_args()
  return args


def main():
  args = parse_args()
  template = 'Reading gammas from {args.gammas}...'
  logging.info(template.format(**locals()))
  data = pd.read_csv(args.gammas, sep='\t')
  logging.info('Drawing plot...')
  fig, ax = plt.subplots(1, 1, figsize=(6,6))
  subs = data[['variant', 'rep', 'gamma']]
  subs = subs.set_index(['variant', 'rep'])
  subs = subs.unstack(level='rep')
  a = subs.columns[-2]
  b = subs.columns[-1]
  plt.xlim(-1.3, 0.1)
  plt.ylim(-1.3, 0.1)
  plot = sns.scatterplot(a, b, data=subs,
                         s=5, alpha=0.5, edgecolor='none')
  plt.xlabel('gamma [replicate a]'.format(**locals()))
  plt.ylabel('gamma [replicate b]'.format(**locals()))
  cleansubs = subs.dropna(subset=[a, b])
  prs, _ = st.pearsonr(cleansubs[a], cleansubs[b])
  ax.text(0.05, 0.95,
          'Pearson R: {prs:.2f}'.format(**locals()),
          transform=ax.transAxes)

  plt.tight_layout()
  plt.savefig(args.pngfile, dpi=600)
  plt.close('all')

##############################################
if __name__ == "__main__":
  sys.exit(main())

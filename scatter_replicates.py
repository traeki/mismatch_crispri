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
      '--organism', type=str,
      help='str: name of the organism',
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
  figure = plt.figure(figsize=(6,6))
  subs = data[['variant', 'rep', 'gamma']]
  subs = subs.set_index(['variant', 'rep'])
  subs = subs.unstack(level='rep')
  a = subs.columns[-2]
  b = subs.columns[-1]
  template = 'Comparison of representative replicates for {args.organism}'
  plt.title(template.format(**vars()))
  plt.xlim(-1.3, 0.1)
  plt.ylim(-1.3, 0.1)
  plot = sns.scatterplot(a, b, data=subs,
                         s=5, alpha=0.5, edgecolor='none')
  plt.xlabel('gammas measured in replicate {a[1]}'.format(**locals()))
  plt.ylabel('gammas measured in replicate {b[1]}'.format(**locals()))
  plt.tight_layout()
  plt.savefig(args.pngfile, dpi=300)
  plt.close('all')

##############################################
if __name__ == "__main__":
  sys.exit(main())

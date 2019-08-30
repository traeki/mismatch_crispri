#!/usr/bin/env python
# Author: John Hawkins (jsh) [really@gmail.com]

import argparse
import logging
import pathlib
import shutil
import sys

from matplotlib import pyplot as plt
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
  prs = el.pearsons_by_gene(data)

  fig = plt.figure(figsize=(8,4))
  sns.distplot(prs.dropna(), norm_hist=True)
  medprs = prs.median()

  plt.text(0, 1, 'median: {medprs}'.format(**locals()))
  template = 'Pearson Correlation Distribution for {args.name}'
  plt.title(template.format(**vars()))
  plt.xlim(-0.1, 1.1)
  plt.xlabel(args.name)
  plt.tight_layout()
  plt.savefig(args.pngfile, dpi=300)
  plt.close('all')

##############################################
if __name__ == "__main__":
  sys.exit(main())

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
      '--xrelfits', type=str,
      help='file: unique entries for relfit by variant -- x axis',
      required=True)
  parser.add_argument(
      '--xname', type=str,
      help='str: label name for x axis',
      required=True)
  parser.add_argument(
      '--yrelfits', type=str,
      help='file: unique entries for relfit by variant -- y axis',
      required=True)
  parser.add_argument(
      '--yname', type=str,
      help='str: label name for y axis',
      required=True)
  parser.add_argument(
      '--pngfile', type=str,
      help='directory: directory for plots (WARNING: will be created and cleared)',
      required=True)
  args = parser.parse_args()
  return args


def main():
  args = parse_args()
  template = 'Reading X: {args.xname} from {args.xrelfits}...'
  logging.info(template.format(**locals()))
  xdata = pd.read_csv(args.xrelfits, sep='\t')
  xprs = el.pearsons_by_gene(xdata)

  template = 'Reading Y: {args.yname} from {args.yrelfits}...'
  logging.info(template.format(**locals()))
  ydata = pd.read_csv(args.yrelfits, sep='\t')
  yprs = el.pearsons_by_gene(ydata)

  hues = ydata.groupby('gene').relfit.min().abs()

  data = pd.concat([xprs, yprs], axis='columns')
  data.columns = [args.xname, args.yname]

  fig = plt.figure(figsize=(6,6))
  sns.scatterplot(args.xname, args.yname, data=data,
                  hue=hues, legend=False,
                  s=10, alpha=1, edgecolor='none')

  plt.xlim(-0.1, 1.1)
  plt.ylim(-0.1, 1.1)
  template = 'Pearson R for each gene in {0}'
  plt.xlabel(template.format(args.xname))
  plt.ylabel(template.format(args.yname))
  plt.tight_layout()
  plt.savefig(args.pngfile, dpi=600)
  plt.close('all')

##############################################
if __name__ == "__main__":
  sys.exit(main())

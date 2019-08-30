#!/usr/bin/env python
# Author: John Hawkins (jsh) [really@gmail.com]

import argparse
import logging
import pandas as pd
import pathlib
import shutil
import sys

from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s') 

_PACKAGEDIR = pathlib.Path(__file__).parent
TESTDIR = _PACKAGEDIR / 'testdata'
_CODEFILE = pathlib.Path(__file__)
NORMCOUNT = 100*1000*1000


def parse_args():
  """Read in the arguments for the sgrna library construction code."""
  logging.info('Parsing command line.')
  parser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument(
      '--xcounts', type=str,
      help='file: unique entries for count by variant -- x axis',
      required=True)
  parser.add_argument(
      '--xname', type=str,
      help='str: label name for x axis',
      required=True)
  parser.add_argument(
      '--ycounts', type=str,
      help='file: unique entries for count by variant -- y axis',
      required=True)
  parser.add_argument(
      '--yname', type=str,
      help='str: label name for y axis/legend',
      required=True)
  parser.add_argument(
      '--guide_subset', type=str,
      help='file: list of variants to include',
      required=False)
  parser.add_argument(
      '--pngfile', type=str,
      help='file: output destination for plot',
      required=True)
  args = parser.parse_args()
  return args


def main():
  args = parse_args()

  template = 'Reading X: {args.xname} from {args.xcounts}...'
  logging.info(template.format(**locals()))
  xdata = pd.read_csv(args.xcounts, sep='\t', header=None)
  xdata.columns = ['variant', args.xname]
  xdata[args.xname] = xdata[args.xname] * NORMCOUNT / xdata[args.xname].sum()
  xdata[args.xname] = xdata[args.xname].apply(np.log2).clip(0)
  xmax = xdata[args.xname].max()

  template = 'Reading Y: {args.yname} from {args.ycounts}...'
  logging.info(template.format(**locals()))
  ydata = pd.read_csv(args.ycounts, sep='\t', header=None)
  ydata.columns = ['variant', args.yname]
  ydata[args.yname] = ydata[args.yname] * NORMCOUNT / ydata[args.yname].sum()
  ydata[args.yname] = ydata[args.yname].apply(np.log2).clip(0)
  ymax = ydata[args.yname].max()

  data = pd.merge(xdata, ydata, on='variant', how='inner')
  if args.guide_subset is not None:
    subset = set(pd.read_csv(args.guide_subset, header=None)[0])
    data = data.loc[data.variant.isin(subset)]
  datamax = max(xmax, ymax)
  data = data.dropna(how='any')
  logging.info('Drawing plot to {args.pngfile}...'.format(**locals()))
  figure = plt.figure(figsize=(6,6))
  plot = sns.scatterplot(args.xname, args.yname, data=data,
                         s=5, alpha=0.5, edgecolor='none')
  # plt.title('Raw Count\n{args.xname}\nvs\n{args.yname}'.format(**vars()))
  plt.xlim(-0.2, datamax)
  plt.ylim(-0.2, datamax)
  plt.xlabel(args.xname)
  plt.ylabel(args.yname)
  plt.tight_layout()
  plt.savefig(args.pngfile, dpi=300)
  plt.close('all')

##############################################
if __name__ == "__main__":
  sys.exit(main())

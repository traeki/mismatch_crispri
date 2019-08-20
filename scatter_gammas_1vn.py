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
  """Read in the arguments for the sgrna library construction code."""
  logging.info('Parsing command line.')
  parser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument(
      '--xgammas', type=str,
      help='file: unique entries for gamma by variant -- x axis',
      required=True)
  parser.add_argument(
      '--xname', type=str,
      help='str: label name for x axis',
      required=True)
  parser.add_argument(
      '--ygammas', type=str,
      action='append',
      help='file: unique entries for gamma by variant -- y axis (repeatable)',
      required=True)
  parser.add_argument(
      '--yname', type=str,
      action='append',
      help='str: label name for y axis/legend',
      required=True)
  parser.add_argument(
      '--pngfile', type=str,
      help='directory: directory for plots (WARNING: will be created and cleared)',
      required=True)
  args = parser.parse_args()
  assert len(args.ygammas) == len(args.yname)
  return args


def main():
  args = parse_args()
  template = 'Reading X: {args.xname} from {args.xgammas}...'
  logging.info(template.format(**locals()))
  xdata = pd.read_csv(args.xgammas, sep='\t')
  data = xdata[['variant', 'gamma', 'gene']]
  data.columns = ['variant', 'base', 'gene']
  hues = list()
  for yname, ygammas in zip(args.yname, args.ygammas):
    template = 'Reading Y: {args.yname} from {args.ygammas}...'
    logging.info(template.format(**locals()))
    yg = pd.read_csv(ygammas, sep='\t')
    yg = yg[['variant', 'gamma']]
    yd = pd.merge(data, yg, on='variant', how='outer')
    yd = yd.dropna(how='any', subset=['base', 'gamma'])
    yd['name'] = yname
    hues.append(yd)
  ydata = pd.concat(hues, axis='index')
  logging.info('Drawing plot...')
  figure = plt.figure(figsize=(6,6))
  plot = sns.scatterplot('base', 'gamma', data=ydata, hue='name',
                         s=5, alpha=0.5, edgecolor='none')
  flatnames = ', '.join(['{0}'.format(x) for x in args.yname])
  plt.title('{args.xname}\nvs\n{flatnames}'.format(**vars()))
  plt.xlim(-1.3, 0.1)
  plt.ylim(-1.3, 0.1)
  plt.xlabel(args.xname)
  plt.ylabel('')
  plt.legend(loc='lower left', title=None)
  plt.tight_layout()
  plt.savefig(args.pngfile, dpi=300)
  plt.close('all')

##############################################
if __name__ == "__main__":
  sys.exit(main())

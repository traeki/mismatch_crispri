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

GAMMARANGE = (-1.3, 0.1)

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
      '--controls', type=str,
      help='file: list of control guides to exclude',
      default=None)
  parser.add_argument(
      '--plotdir', type=str,
      help='directory: directory for plots (WARNING: will be created and cleared)',
      required=True)
  parser.add_argument(
      '--pngfile', type=str,
      help='str: template for filenames to write inside plotdir',
      default='1vn.scatter.png')
  args = parser.parse_args()
  if len(args.ygammas) != len(args.yname):
    logging.error('Must provide --yname for every instance of --ygammas.')
    sys.exit(3)
  return args


def main():
  args = parse_args()
  if args.controls is not None:
    controls = set(pd.read_csv(args.controls, header=None)[0])
  else:
    controls = set()
  template = 'Reading X: {args.xname} from {args.xgammas}...'
  logging.info(template.format(**locals()))
  xdata = pd.read_csv(args.xgammas, sep='\t')
  data = xdata[['variant', 'gamma', 'gene']]
  data.columns = ['variant', 'base', 'gene']
  hues = list()
  for yname, ygammas in zip(args.yname, args.ygammas):
    template = 'Reading Y: {yname} from {ygammas}...'
    logging.info(template.format(**locals()))
    yg = pd.read_csv(ygammas, sep='\t')
    yg = yg[['variant', 'gamma']]
    yd = pd.merge(data, yg, on='variant', how='outer')
    yd = yd.dropna(how='any', subset=['base', 'gamma'])
    yd['name'] = yname
    hues.append(yd)
  ydata = pd.concat(hues, axis='index')
  ydata = ydata.loc[~ydata.variant.isin(controls)]

  plotdir = pathlib.Path(args.plotdir)
  shutil.rmtree(plotdir, ignore_errors=True)
  plotdir.mkdir(parents=True, exist_ok=True)
  for gene, group in ydata.groupby('gene'):
    suffix = '.{gene}.png'.format(**locals())
    plotfile = (plotdir / args.pngfile).with_suffix(suffix)
    logging.info('Drawing plot for {gene}...'.format(**locals()))
    figure = plt.figure(figsize=(6,6))
    plt.xlim(*GAMMARANGE)
    plt.ylim(*GAMMARANGE)
    if len(args.yname) == 1:
      plot = sns.scatterplot('base', 'gamma', data=group, hue='name',
                             s=10, alpha=0.5, edgecolor='none', legend=False)
      plt.ylabel(args.yname[0])
    else:
      plot = sns.scatterplot('base', 'gamma', data=group, hue='name',
                             s=10, alpha=0.5, edgecolor='none', legend='brief')
      plt.ylabel('')
    plt.plot(GAMMARANGE, GAMMARANGE, 'b--', linewidth=0.5)
    plt.xlabel(args.xname)
    plt.title(gene)
    plt.tight_layout()
    plt.savefig(plotfile, dpi=600)
    plt.close('all')

##############################################
if __name__ == "__main__":
  sys.exit(main())

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

RELFITRANGE = (-0.3, 1.1)

def parse_args():
  """Read in the arguments for the sgrna library construction code."""
  logging.info('Parsing command line.')
  parser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument(
      '--xrelfits', type=str,
      help='file: unique entries for fitness by variant -- x axis',
      required=True)
  parser.add_argument(
      '--xname', type=str,
      help='str: label name for x axis',
      required=True)
  parser.add_argument(
      '--yrelfits', type=str,
      action='append',
      help='file: unique entries for fitness by variant -- y axis (repeatable)',
      required=True)
  parser.add_argument(
      '--yname', type=str,
      action='append',
      help='str: label name for y axis/legend',
      required=True)
  parser.add_argument(
      '--ylabel', type=str,
      help='str: override label for y axis',
      default=None)
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
  if len(args.yrelfits) != len(args.yname):
    logging.error('Must provide --yname for every instance of --yrelfits.')
    sys.exit(3)
  return args


def main():
  args = parse_args()
  if args.controls is not None:
    controls = set(pd.read_csv(args.controls, header=None)[0])
  else:
    controls = set()
  template = 'Reading X: {args.xname} from {args.xrelfits}...'
  logging.info(template.format(**locals()))
  xdata = pd.read_csv(args.xrelfits, sep='\t')
  data = xdata[['variant', 'relfit', 'gene']]
  data.columns = ['variant', 'base', 'gene']
  hues = list()
  for yname, yrelfits in zip(args.yname, args.yrelfits):
    template = 'Reading Y: {yname} from {yrelfits}...'
    logging.info(template.format(**locals()))
    yg = pd.read_csv(yrelfits, sep='\t')
    yg = yg[['variant', 'relfit']]
    yd = pd.merge(data, yg, on='variant', how='outer')
    yd = yd.dropna(how='any', subset=['base', 'relfit'])
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
    plt.xlim(*RELFITRANGE)
    plt.ylim(*RELFITRANGE)
    if args.ylabel is not None:
      plot = sns.scatterplot('base', 'relfit', data=group, hue='name',
                             s=10, alpha=0.5, edgecolor='none', legend='brief')
      plt.ylabel(args.ylabel)
    elif len(args.yname) == 1:
      plot = sns.scatterplot('base', 'relfit', data=group, hue='name',
                             s=10, alpha=0.5, edgecolor='none', legend=False)
      plt.ylabel(args.yname[0])
    else:
      plot = sns.scatterplot('base', 'relfit', data=group, hue='name',
                             s=10, alpha=0.5, edgecolor='none', legend='brief')
      plt.ylabel('')
    plt.plot(RELFITRANGE, RELFITRANGE, 'b--', linewidth=0.5)
    plt.xlabel(args.xname)
    plt.title(gene)
    plt.tight_layout()
    plt.savefig(plotfile, dpi=600)
    plt.close('all')

##############################################
if __name__ == "__main__":
  sys.exit(main())

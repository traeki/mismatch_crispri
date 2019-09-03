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
      '--basegammas', type=str,
      help='file: unique entries for gamma by variant -- base',
      required=True)
  parser.add_argument(
      '--basename', type=str,
      help='str: label name for base set',
      required=True)
  parser.add_argument(
      '--compgammas', type=str,
      action='append',
      help='file: unique entries for gamma by variant -- comparison (repeatable)',
      required=True)
  parser.add_argument(
      '--compname', type=str,
      action='append',
      help='str: label name for comparison set',
      required=True)
  parser.add_argument(
      '--controls', type=str,
      help='file: list of control guides',
      required=True)
  parser.add_argument(
      '--pngfile', type=str,
      help='file: output location for violin plot',
      default='fig_kernel_violin.png')
  args = parser.parse_args()
  return args


def main():
  args = parse_args()
  if args.controls is not None:
    controls = set(pd.read_csv(args.controls, header=None)[0])
  else:
    controls = set()
  template = 'Reading base set: {args.basename} from {args.basegammas}...'
  logging.info(template.format(**locals()))
  basedata = pd.read_csv(args.basegammas, sep='\t')
  data = basedata[['variant', 'gamma', 'gene']]
  data.columns = ['variant', 'base', 'gene']
  hues = list()
  for compname, compgammas in zip(args.compname, args.compgammas):
    template = 'Reading comp: {compname} from {compgammas}...'
    logging.info(template.format(**locals()))
    compg = pd.read_csv(compgammas, sep='\t')
    compg = compg[['variant', 'gamma']]
    compd = pd.merge(data, compg, on='variant', how='outer')
    compd = compd.dropna(how='any', subset=['base', 'gamma'])
    compd['name'] = compname
    hues.append(compd)
  compdata = pd.concat(hues, axis='index')
  cmask = compdata.variant.isin(controls)
  controldata = compdata.loc[cmask]
  compdata = compdata.loc[~cmask]
  radius = controldata.base.std()
  nullset = compdata.loc[compdata.base.abs() < radius]

  figure = plt.figure(figsize=(6,6))
  logging.info('Drawing plot to {args.pngfile}...'.format(**locals()))

  subset = nullset.loc[nullset.gene == 'dfrA']
  sns.violinplot(x='name', y='gamma', data=subset, scale='width')

  plt.tight_layout()
  plt.savefig(args.pngfile, dpi=600)
  plt.close('all')

##############################################
if __name__ == "__main__":
  sys.exit(main())

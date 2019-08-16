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


def parse_args():
  """Read in the arguments for the sgrna library construction code."""
  logging.info('Parsing command line.')
  parser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument(
      '--meangammas', type=str,
      help='file: file containing averaged gammas (i.e. ???.gammas.mean.tsv)',
      default=str(TESTDIR / 'test.config.tsv.gammas.mean.tsv'))
  parser.add_argument(
      '--plotdir', type=str,
      help='directory: directory for plots (WARNING: will be created and cleared)',
      default=str(TESTDIR / 'gvkplots'))
  args = parser.parse_args()
  return args


def main():
  args = parse_args()
  # reset PLOTDIR
  data = pd.read_csv(args.meangammas, sep='\t')
  data.set_index('variant', inplace=True)
  data['knockdown'] = data['y_pred']
  plotdir = pathlib.Path(args.plotdir)
  shutil.rmtree(plotdir, ignore_errors=True)
  plotdir.mkdir(parents=True, exist_ok=True)
  # draw gene-by-gene scatterplots
  logging.info('Drawing plots...')
  for gene, group in data.groupby('gene'):
    plotfile = plotdir / '.'.join(['gvk', gene, 'png'])
    figure = plt.figure(figsize=(6,6))
    plot = sns.scatterplot('knockdown', 'gamma', data=group, hue='original',
                           s=10, alpha=1, edgecolor='none', legend=False)
    plt.title('{gene}\nKnockdown vs. Gamma'.format(**vars()))
    plt.xlim(-0.1, 1.1)
    plt.ylim(-1.3, 0.1)
    plt.xlabel('Knockdown (predicted)')
    plt.ylabel('Pooled-growth gamma')
    plt.tight_layout()
    plt.savefig(plotfile, dpi=300)
    plt.close('all')

##############################################
if __name__ == "__main__":
  sys.exit(main())

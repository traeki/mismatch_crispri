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
import scipy.stats as st
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
      '--meanrelfit', type=str,
      help='file: file containing averaged relfit (e.g. relfit.mean.tsv)',
      default=str(TESTDIR / 'relfit.mean.tsv'))
  parser.add_argument(
      '--plotdir', type=str,
      help='directory: directory for plots (WARNING: will be created and cleared)',
      default=str(TESTDIR / 'kvf.plots'))
  args = parser.parse_args()
  return args


def plot_kvf(data, name, plotfile, *, color=True):
  data = data.dropna(subset=['knockdown', 'relfit'])
  if len(data) < 1:
    logging.warn('No data to plot for {name}'.format(**locals()))
    return
  if len(data) > 1:
    prs, _ = st.pearsonr(data.knockdown, data.relfit)
  else:
    prs = np.nan
  figure = plt.figure(figsize=(6,6))
  hue = (color and 'original' or None)
  plot = sns.scatterplot('knockdown', 'relfit', data=data, hue=hue,
                         s=10, alpha=1, edgecolor='none', legend=False)
  plt.text(0, -1.1, 'Pearson R: {prs:.2f}'.format(**locals()))
  plt.title('{name}\nKnockdown vs. Relative Fitness'.format(**vars()))
>>>>>>> master:kvf_by_gene.py
  plt.xlim(-0.1, 1.1)
  plt.ylim(-0.3, 1.1)
  plt.xlabel('Knockdown (predicted)')
  plt.ylabel('Relative Pooled-growth Fitness')
  plt.tight_layout()
  plt.savefig(plotfile, dpi=600)
  plt.close('all')

def main():
  args = parse_args()
  # reset PLOTDIR
  data = pd.read_csv(args.meanrelfit, sep='\t')
  data.set_index('variant', inplace=True)
  data['knockdown'] = data['y_pred']
  plotdir = pathlib.Path(args.plotdir)
  shutil.rmtree(plotdir, ignore_errors=True)
  plotdir.mkdir(parents=True, exist_ok=True)
  # draw gene-by-gene scatterplots
  logging.info('Drawing plots...')
  plotfile = plotdir / '.'.join(['kvf', 'overall', 'png'])
  plot_kvf(data, 'OVERALL', plotfile, color=False)
  for gene, group in data.groupby('gene'):
    plotfile = plotdir / '.'.join(['kvf', gene, 'png'])
    plot_kvf(group, gene, plotfile)

##############################################
if __name__ == "__main__":
  sys.exit(main())

#!/usr/bin/env python
# Author: John Hawkins (jsh) [really@gmail.com]

import logging
import numpy as np
import pandas as pd
import pathlib
import shutil

from matplotlib import pyplot as plt
import seaborn as sns

import model_lib as ml

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s') 

_CODEFILE = pathlib.Path(__file__).name
UNGD = pathlib.Path('/home/jsh/ungd/proj/smidge')
DATADIR = UNGD / 'bsubox'
GAMMAFILE = DATADIR / 'compute_gammas.bsu.explore.tsv'
PLOTDIR = (DATADIR / _CODEFILE).with_suffix('.plots')

data = pd.read_csv(GAMMAFILE, sep='\t')
logging.info('Predicting outcomes...')
data['y_pred'] = -(ml.predict_mismatch_scores(data))
anno = data.drop(['rep', 'gamma', 'start_mask'], axis='columns')
anno = anno.drop_duplicates()
gamma = data[['variant', 'rep', 'gamma']].set_index(['variant', 'rep'])
gamma = gamma.unstack(level='rep').mean(axis='columns')
data = pd.DataFrame(anno)
data.set_index('variant', inplace=True)
data['gamma'] = gamma

# reset PLOTDIR
shutil.rmtree(PLOTDIR, ignore_errors=True)
PLOTDIR.mkdir(parents=True, exist_ok=True)
# draw gene-by-gene scatterplots
logging.info('Drawing plots...')
for gene, group in data.groupby('gene'):
  plotfile = PLOTDIR / '.'.join(['pvm', gene, 'png'])
  figure = plt.figure(figsize=(6,6))
  plot = sns.scatterplot('y_pred', 'gamma', data=group, hue='original',
                         s=10, alpha=1, edgecolor='none', legend=False)
  plt.title('{gene}\nPredicted vs. Measured'.format(**vars()))
  plt.xlim(-1.1, 0.1)
  plt.ylim(-1.1, 0.1)
  plt.tight_layout()
  plt.savefig(plotfile, dpi=300)
  plt.close('all')

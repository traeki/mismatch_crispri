#!/usr/bin/env python
# Author: John Hawkins (jsh) [really@gmail.com]

import logging
import sys

import pandas as pd
import numpy as np
import scipy.stats as st
import seaborn as sns

import choice_lib as cl

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

def pearsons_by_gene(frame):
  pearsons = dict()
  for gene, group in frame.groupby('gene'):
    good = group.dropna(how='any', subset=['y_pred', 'gamma'])
    predicted = good.y_pred
    measured = -good.gamma
    try:
      prsrho, prspv = st.pearsonr(predicted, measured)
      pearsons[gene] = prsrho
    except ValueError:
      pearsons[gene] = np.nan
  result = pd.Series(pearsons)
  return result

def bin_counts_by_gene(frame):
  bin_edges = cl.gamma_bins()
  frame['bin'] = cl.bin_gammas(frame.gamma, bin_edges)
  bintally = frame.groupby('gene').bin.value_counts().unstack()
  bintally['mid'] = bintally[[1, 2, 3]].sum(axis='columns')
  return bintally

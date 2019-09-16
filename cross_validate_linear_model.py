#!/usr/bin/env python
# Author: John Hawkins (jsh) [really@gmail.com]

import joblib
import logging
import pathlib
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as st

import model_lib as ml

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

_CODEDIR = pathlib.Path(__file__).parent
MODELDIR = _CODEDIR / 'model'
GFPDIR = _CODEDIR / 'gfpdata'

bsugfp = pd.read_csv(GFPDIR / 'gfp.bsu.facsseq.tsv', sep='\t')
ecogfp = pd.read_csv(GFPDIR / 'gfp.eco.facsseq.tsv', sep='\t')
origmap = pd.read_csv(GFPDIR / 'gfp.origmap.tsv', sep='\t')

nmm = origmap[['nmm']]
origmap = origmap[['variant', 'original']]

def relscore(gfpscore, vomap):
  gfpscore = gfpscore.set_index('variant').facsseq
  def lookup(x):
    if x in gfpscore:
      return gfpscore.loc[x]
    else:
      return np.nan
  scoremap = vomap.applymap(lookup)
  return (scoremap.variant / scoremap.original)

data = pd.DataFrame(origmap)
bsu_rs = relscore(bsugfp, origmap)
eco_rs = relscore(ecogfp, origmap)
score = pd.concat([bsu_rs, eco_rs], axis='columns').mean(axis='columns')
data['y'] = score
data = data.loc[nmm.nmm == 1]
data = data.dropna(axis='rows')
mm_data = data[['variant', 'original']]
y_data = data[['y']]

ml.kfold_score_model(mm_data, y_data)

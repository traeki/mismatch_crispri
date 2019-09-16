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

repfiles = ['bsu_biorep1.csv',
            'bsu_biorep2.csv',
            'eco_biorep1.csv',
            'eco_biorep2.csv']
replicates = list()
for repfile in repfiles:
  repdata = pd.read_csv(GFPDIR / repfile, index_col=0)
  repdata = repdata[['relative']].dropna()
  replicates.append(repdata)
score = pd.concat(replicates, axis='columns', sort=True).mean(axis='columns')
# "relative" in these files is (C/P)-1 ; downstream assumes C/P
score = score + 1
score = pd.DataFrame(score).reset_index()
score.columns = ['variant', 'y']

origmap = pd.read_csv(GFPDIR / 'gfp.origmap.tsv', sep='\t')
nmm = origmap[['nmm']]
origmap = origmap[['variant', 'original']]
data = pd.DataFrame(origmap)
data = data.merge(score, on='variant', how='left')

data = data.loc[nmm.nmm == 1]
data = data.dropna(axis='rows')
mm_data = data[['variant', 'original']]
y_data = data[['y']]

ml.train_and_save_mismatch_model(mm_data, y_data)

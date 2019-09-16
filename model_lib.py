#!/usr/bin/env python
# Author: John Hawkins (jsh) [really@gmail.com]

import itertools
import joblib
import logging
import os.path
import pathlib
import random
import shutil
import sys

import numpy as np
import pandas as pd
from pandas.api.types import CategoricalDtype
from sklearn import preprocessing as skpreproc
from keras.layers import Dense
from keras.models import Sequential


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')
np.set_printoptions(precision=4, suppress=True)


_CODEDIR = pathlib.Path(__file__).parent
MODELDIR = _CODEDIR / 'model'

MODELFILE = MODELDIR / 'model.d5'.format(**locals())
XS_FILE = MODELDIR / 'xscaler.dump'
YS_FILE = MODELDIR / 'yscaler.dump'

_EPOCHS = 30
_BATCH_SIZE = 32

def _build_linear_model(num_features):
  model = Sequential()
  model.add(Dense(1, input_dim=num_features, activation='linear'))
  model.compile(loss='mse', metrics=['mse'], optimizer='adam')
  return model

def _expand_dummies(frame):
  categories = dict()
  bases = ['A', 'C', 'G', 'T']
  idxs = [x for x in range(20)]  # Magic number because guidelen is fixed.
  pairs = [''.join(pair) for pair in itertools.product(bases, bases)]
  categories['mm_idx'] = idxs
  categories['mm_trans'] = pairs
  widecols = list()
  for column in frame.columns:
    if column not in categories:
      continue
    frame[column] = frame[column].astype(CategoricalDtype(categories[column]))
  return pd.get_dummies(frame)

def _get_linear_encoder():
  def encoder(inrow):
    vari = inrow.variant
    orig = inrow.original
    mm_idx = None
    for i in range(len(vari)):
      if vari[i] != orig[i]:
        if orig != vari[:i] + orig[i] + vari[i+1:]:
          template = 'too many mismatches in pair {vari} <- {orig}'
          raise ValueError(template.format(**locals()))
        mm_idx = i
    if mm_idx == None:
      template = 'no mismatch in pair {vari} <- {orig}'
      raise ValueError(template.format(**locals()))
    features = dict()
    features['mm_idx'] = mm_idx
    mm_trans = ''.join([orig[mm_idx], vari[mm_idx]])
    features['mm_trans'] = mm_trans
    features['gc_cont'] = orig.count('G') + orig.count('C')
    row = pd.Series(features)
    return row
  return encoder

def train_and_save_mismatch_model(voframe, yframe):
  if voframe.shape[0] != yframe.shape[0]:
    logging.fatal('voframe and training values had different length')
    sys.exit(2)
  if 'variant' not in voframe.columns or 'original' not in voframe.columns:
    logging.fatal('voframe missing variant and/or original')
    sys.exit(2)
  encoder = _get_linear_encoder()
  encodings = voframe.apply(encoder, axis=1)
  Xframe = encodings.set_index(voframe.variant)
  Xframe = _expand_dummies(Xframe)
  X = np.array(Xframe, dtype=float)
  y = np.array(yframe.y, dtype=float).reshape(-1, 1)
  shutil.rmtree(MODELDIR, ignore_errors=True)
  while os.path.exists(MODELDIR):
    continue
  MODELDIR.mkdir(parents=True, exist_ok=True)
  y_orig = y
  X_scaler = skpreproc.StandardScaler()
  X = X_scaler.fit_transform(X)
  y_scaler = skpreproc.StandardScaler()
  y = y_scaler.fit_transform(y)
  model = _build_linear_model(X.shape[1])
  # Feed training Data
  model.fit(X, y, epochs=_EPOCHS, batch_size=_BATCH_SIZE)
  joblib.dump(X_scaler, XS_FILE)
  joblib.dump(y_scaler, YS_FILE)
  joblib.dump(model, MODELFILE)
def _retrieve_mismatch_model():
  try:
    return (joblib.load(MODELFILE), joblib.load(XS_FILE), joblib.load(YS_FILE))
  except FileNotFoundError:
    logging.fatal('Tried to make predictions without a model in place')
    sys.exit(2)

def predict_mismatch_scores(reference):
  if 'variant' not in reference.columns or 'original' not in reference.columns:
    logging.fatal('reference missing variant and/or original')
    sys.exit(2)
  refsize = len(reference)
  logging.info('Applying model to {refsize} guides...'.format(**locals()))
  model, xscaler, yscaler = _retrieve_mismatch_model()
  encoder = _get_linear_encoder()
  voframe = reference[['variant', 'original']]
  voframe = voframe.drop_duplicates()
  matchmask = (voframe.variant == voframe.original)
  parents = pd.DataFrame(voframe.loc[matchmask])
  parents['score'] = 1.0
  children = pd.DataFrame(voframe.loc[~matchmask])
  encodings = children.apply(encoder, axis=1)
  Xframe = encodings.set_index(children.variant)
  Xframe = _expand_dummies(Xframe)
  X = np.array(Xframe, dtype=float)
  X = xscaler.transform(X)
  children['score'] = yscaler.inverse_transform(model.predict(X))
  both = pd.concat([parents, children], axis='rows')
  reconcile = pd.merge(reference, both, on='variant', how='left')
  return reconcile.score

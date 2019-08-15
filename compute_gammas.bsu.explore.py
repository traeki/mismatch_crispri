#!/usr/bin/env python
# Author: John Hawkins (jsh) [really@gmail.com]

import logging
import numpy as np
import pandas as pd
import pathlib
import gamma_lib as gl


_CODEFILE = pathlib.Path(__file__).name
UNGD = pathlib.Path('/home/jsh/ungd/proj/smidge')
DATADIR = UNGD / 'bsubox'
GAMMAFILE = (DATADIR / _CODEFILE).with_suffix('.tsv')

GENBANK = DATADIR / 'bsu.NC_000964.merged.gb'
TARGETFILE = DATADIR / 'bsu.NC_000964.merged.gb.targets.all.tsv'
CONTROLFILE = DATADIR / '2019.chip2.control.guides'
PREDICTFILE = DATADIR / 'predict_all_linear.bsu.tsv'

controlset = gl.get_controlset(CONTROLFILE)
replicates = dict()
replicates['a1'] = (
  DATADIR / 'bsu.explore.A.1.0.counts',
  DATADIR / 'bsu.explore.A.1.10.counts')
replicates['a2'] = (
  DATADIR / 'bsu.explore.A.2.0.counts',
  DATADIR / 'bsu.explore.A.2.10.counts')
replicates['a3'] = (
  DATADIR / 'bsu.explore.A.3.0.counts',
  DATADIR / 'bsu.explore.A.3.10.counts')
replicates['b1'] = (
  DATADIR / 'bsu.explore.B.1.0.counts',
  DATADIR / 'bsu.explore.B.1.10.counts')
replicates['b2'] = (
  DATADIR / 'bsu.explore.B.2.0.counts',
  DATADIR / 'bsu.explore.B.2.10.counts')
replicates['b3'] = (
  DATADIR / 'bsu.explore.B.3.0.counts',
  DATADIR / 'bsu.explore.B.3.10.counts')

rep_frames = list()
for name, (startfile, endfile) in replicates.items():
  gt = 10
  rep_frame = gl.compute_gamma(startfile, endfile, controlset, gt)
  rep_frame['rep'] = name
  rep_frame = rep_frame.reset_index()
  rep_frames.append(rep_frame)
collected = pd.concat(rep_frames).reset_index(drop=True)
variants = collected.variant
annoframe = gl.annotate_variants(variants, PREDICTFILE, TARGETFILE, GENBANK)
annoframe = annoframe.drop('variant', axis='columns')
annoframe = pd.concat([collected, annoframe], axis='columns')
annoframe.to_csv(GAMMAFILE, index=False, sep='\t')

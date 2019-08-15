#!/usr/bin/env python
# Author: John Hawkins (jsh) [really@gmail.com]

import logging
import sys

import pandas as pd
import numpy as np
from Bio import SeqIO

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')

NORM_SIZE = float(100 * 1000 * 1000)
MIN_START_READS = 100
PSEUDO = 1

def get_start_mask(startfile, endfile):
  start = pd.read_csv(startfile, sep='\t', header=None,
                      names=['variant','reads'], index_col='variant')
  end = pd.read_csv(endfile, sep='\t', header=None,
                    names=['variant','reads'], index_col='variant')
  start_mask = start.reads > MIN_START_READS
  start_mask.name = 'start_mask'
  return start_mask

def log_counts(filename):
  reads = pd.read_csv(filename, sep='\t', header=None,
                      names=['variant','reads'], index_col='variant')
  norm = reads.reads * (NORM_SIZE / reads.reads.sum())
  log = np.log2(norm.clip(PSEUDO))
  return log

def get_controlset(controlfile):
  controlframe = pd.read_csv(controlfile, header=None, names=['variant'])
  return set(controlframe.variant)

def compute_gamma(startfile, endfile, controlset, gt):
  start_mask = get_start_mask(startfile, endfile)
  start = log_counts(startfile)
  end = log_counts(endfile)
  diff = end - start
  center = diff.loc[diff.index.isin(controlset)].median()
  gamma = (diff - center) / gt
  gamma.name = 'gamma'
  masked_gamma = gamma.where(start_mask)
  columns = [masked_gamma, start_mask]
  frame = pd.concat(columns, sort=True, axis='columns')
  frame.index.name = 'variant'
  return frame

# TODO(jsh): simplify annotations down to "chosen" file
# TODO(jsh): predictions should be the exlusive purview of model_lib!
def annotate_variants(variants, predictfile, targetfile, genbank):
  annoframe = pd.DataFrame(index=variants)
  predictframe = pd.read_csv(predictfile, sep='\t')
  relevant = predictframe.loc[predictframe.variant.isin(variants) &
                              predictframe.original.isin(variants)]
  relevant = relevant[['variant', 'original']]
  annoframe = pd.merge(annoframe, relevant,
                       left_index=True, right_on='variant', how='left')
  unset_mask = annoframe.original.isna()
  annoframe.original = annoframe.variant.where(unset_mask, annoframe.original)
  origlocusmap = predictframe[['original', 'locus_tag']].drop_duplicates()
  annoframe = pd.merge(annoframe, origlocusmap, on='original', how='left')
  targetframe = pd.read_csv(targetfile, sep='\t', low_memory=False)
  targetframe = targetframe.iloc[:,:4]
  targetframe.columns = ['locus_tag', 'offset', 'original', 'pam']
  annoframe = pd.merge(annoframe, targetframe,
                       on=['original', 'locus_tag'], how='left')
  locusgenemap = dict()
  gbhandle = open(genbank, 'r')
  for record in SeqIO.parse(gbhandle, 'genbank'):
    for feature in record.features:
      if feature.type == 'gene':
        quals = feature.qualifiers
        if 'gene' in quals:
          locusgenemap[quals['locus_tag'][0]] = quals['gene'][0]
        else:
          locusgenemap[quals['locus_tag'][0]] = quals['locus_tag'][0]
  locusgenemap = pd.DataFrame.from_dict(locusgenemap, orient='index')
  locusgenemap = locusgenemap.reset_index()
  locusgenemap.columns = ['locus_tag', 'gene']
  annoframe = pd.merge(annoframe, locusgenemap, on='locus_tag', how='left')
  return annoframe

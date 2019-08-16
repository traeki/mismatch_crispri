#!/usr/bin/env python
# Author: John Hawkins (jsh) [really@gmail.com]

import logging
import sys

import pandas as pd
import numpy as np
from Bio import SeqIO

import choice_lib as cl

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
  diff = diff.where(start_mask, np.nan)
  center = diff.loc[diff.index.isin(controlset)].median()
  gamma = (diff - center) / gt
  gamma.name = 'gamma'
  masked_gamma = gamma.where(start_mask)
  columns = [masked_gamma, start_mask]
  frame = pd.concat(columns, sort=True, axis='columns')
  frame.index.name = 'variant'
  return frame

def annotate_variants(variants, targetfile, locifile, genbank):
  annoframe = pd.DataFrame(index=variants)
  # TODO(jsh): Reading this in is now redundant.  Ugh.  Fix later.
  targetframe = pd.read_csv(targetfile, sep='\t')
  # Still crashing on mismatch shit
  loci = set(pd.read_csv(locifile, sep='\t', header=None)[0])
  variantspace = cl.build_and_filter_pairs(targetfile, loci)
  relevant = variantspace.loc[variantspace.variant.isin(variants) &
                              variantspace.original.isin(variants)]
  relevant = relevant[['variant', 'original']]
  annoframe = pd.merge(annoframe, relevant,
                       left_index=True, right_on='variant', how='left')
  unset_mask = annoframe.original.isna()
  annoframe.original = annoframe.variant.where(unset_mask, annoframe.original)
  origlocusmap = targetframe[['target', 'locus_tag']].drop_duplicates()
  origlocusmap.columns = ['original', 'locus_tag']
  origlocusmap = origlocusmap.loc[origlocusmap.locus_tag.isin(loci)]
  # TODO(jsh): the next line expands annoframe by 1440
  # TODO(jsh): it's not the controls -- 6000 have no locus_tag
  annoframe = pd.merge(annoframe, origlocusmap, on='original', how='left')
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

#!/usr/bin/env python
# Author: John Hawkins (jsh) [really@gmail.com]

import argparse
import logging
import pathlib
import sys

import pandas as pd

import gamma_lib as gl
from matplotlib import pyplot as plt
import model_lib as ml
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
      '--scorefile', type=str,
      help='file: tsv file of scored variants (flag repeats are averaged)',
      action='append',
      required=True)
  parser.add_argument(
      '--origmap', type=str,
      help='file: variant-original map to override single-mismatch filter',
      required=True)
  parser.add_argument(
      '--pngfile', type=str,
      help='file: output file',
      required=True)
  args = parser.parse_args()
  return args

def count_mismatches(row):
  n = 0
  v = row.variant
  o = row.original
  for i in range(len(v)):
    if v[i] != o[i]:
      n += 1
  return n

def is_double_mismatch(row):
  return count_mismatches(row) == 2

def double_mismatches(frame):
  return frame.apply(is_double_mismatch, axis='columns')

def sub_variants(row):
  o = row.original
  v = row.variant
  subs = list()
  for i in range(len(v)):
    if v[i] != o[i]:
      s = o[:i] + v[i] + o[i+1:]
      subs.append(s)
  assert len(subs) == 2
  return pd.Series({'vara':subs[0], 'varb':subs[1]})

def bubble_gap(row):
  o = row.original
  v = row.variant
  locs = list()
  for i in range(len(v)):
    if v[i] != o[i]:
      locs.append(i)
  assert len(locs) == 2
  return locs[1] - locs[0] - 1

def plot_gvk_doubles(data, compcol, pngfile):
  figure = plt.figure(figsize=(6,6))
  data = data.dropna(subset=[compcol, 'score'])
  prs, pval = st.pearsonr(data[compcol], data.score)
  plot = sns.scatterplot(compcol, 'score', data=data, hue='gap',
                         s=10, alpha=1, edgecolor='none', legend=False)
  plt.text(0, 1.0, 'Pearson R: {prs:.2f}'.format(**locals()))
  plt.xlim(-0.1, 1.1)
  plt.ylim(-0.1, 1.1)
  template = 'Knockdown (synthesis of predictions -- {compcol})'
  plt.xlabel(template.format(**locals()))
  plt.ylabel('FACS-seq score')
  plt.tight_layout()
  logging.info('Drawing {compcol} eval to {pngfile}...'.format(**locals()))
  plt.savefig(pngfile, dpi=600)
  plt.close('all')

def main():
  args = parse_args()
  scoresets = list()
  for scorefile in args.scorefile:
    ss = pd.read_csv(scorefile, sep='\t')
    ss.columns = ['variant', 'gamma']
    ss = ss.set_index('variant')
    scoresets.append(ss)
  flatframe = pd.concat(scoresets, axis='columns', sort=True)
  flatframe = flatframe.mean(axis='columns')
  flatframe = flatframe.reset_index()
  flatframe.columns = ['variant', 'gamma']
  origmap = pd.read_csv(args.origmap, sep='\t').drop_duplicates()
  flatframe = flatframe.merge(origmap, on='variant')
  flatframe = flatframe.loc[double_mismatches(flatframe)]
  subvars = flatframe.apply(sub_variants, axis='columns')
  flatframe = pd.concat([flatframe, subvars], axis='columns')
  flatframe = flatframe.reset_index(drop=True)
  aref = pd.DataFrame(flatframe[['original', 'vara']])
  aref.columns = ['original', 'variant']
  flatframe['a_pred'] = ml.predict_mismatch_scores(aref)
  bref = pd.DataFrame(flatframe[['original', 'varb']])
  bref.columns = ['original', 'variant']
  flatframe['b_pred'] = ml.predict_mismatch_scores(bref)
  flatframe['geometric'] = flatframe.a_pred * flatframe.b_pred
  flatframe['gap'] = flatframe.apply(bubble_gap, axis='columns')
  flatframe['score'] = (-1 * flatframe.gamma)
  compcol = 'geometric'
  suffix = '.{compcol}.png'.format(**locals())
  pngfile = pathlib.Path(args.pngfile).with_suffix(suffix)
  plot_gvk_doubles(flatframe, compcol, pngfile)

##############################################
if __name__ == "__main__":
  sys.exit(main())

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
      '--genbank', type=str,
      help='file: base genome for organism in genbank format',
      default=str(TESTDIR / 'bsu.NC_000964.merged.gb'))
  parser.add_argument(
      '--targetfile', type=str,
      help='file: ???.targets.all.tsv produced by traeki/sgrna_design',
      default=str(TESTDIR / 'test.gb.targets.all.tsv'))
  parser.add_argument(
      '--controls', type=str,
      help='file: list of just control guides',
      default=str(TESTDIR / 'test.controls'))
  parser.add_argument(
      '--locifile', type=str,
      help='file: list of applicable locus_tags',
      default=str(TESTDIR / 'test.loci'))
  parser.add_argument(
      '--justgene', type=str,
      help='str: gene to graph in isolation (default is to plot all genes)',
      default=None)
  parser.add_argument(
      '--configdir', type=str, help='file: name of directory containing config.tsv',
      default=str(TESTDIR))
  parser.add_argument(
      '--origmap', type=str,
      help='file: variant-original map to override single-mismatch filter',
      required=True)
  parser.add_argument(
      '--pngfile', type=str,
      help='file: output file',
      required=True)
  parser.add_argument(
      '--growth', type=int,
      help='int: number of generations grown (in other words, g*t)',
      default=10)
  args = parser.parse_args()
  return args


def flatgamma(stacked_replicates, controls):
  data = stacked_replicates
  anno = data.drop(['rep', 'gamma', 'start_mask'], axis='columns')
  anno = anno.drop_duplicates()
  gamma = data[['variant', 'rep', 'gamma']].set_index(['variant', 'rep'])
  gamma = gamma.unstack(level='rep').mean(axis='columns')
  data = pd.DataFrame(anno)
  data.set_index('variant', inplace=True)
  data['gamma'] = gamma
  return data


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

def n_mismatches(frame):
  return frame.apply(count_mismatches, axis='columns')

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
  data = data.dropna(subset=[compcol, 'gamma'])
  prs, pval = st.pearsonr(data[compcol], -data.gamma)
  plot = sns.scatterplot(compcol, 'gamma', data=data, hue='gap',
                         s=10, alpha=1, edgecolor='none', legend=False)
  plt.text(0, -1.1, 'Pearson R: {prs:.2f}'.format(**locals()))
  plt.xlim(-0.1, 1.1)
  plt.ylim(-1.3, 0.1)
  template = 'Knockdown (synthesis of predictions -- {compcol})'
  plt.xlabel(template.format(**locals()))
  plt.ylabel('Pooled-growth gamma')
  plt.tight_layout()
  logging.info('Drawing {compcol} eval to {pngfile}...'.format(**locals()))
  plt.savefig(pngfile, dpi=600)
  plt.close('all')

def main():
  args = parse_args()
  controls = set(pd.read_csv(args.controls, header=None)[0])
  configdir = pathlib.Path(args.configdir)
  config = pd.read_csv(configdir / 'config.tsv', sep='\t')
  config = config.set_index('sample')
  rep_frames = list()
  for rep, ends in config.iterrows():
    logging.info('Computing gammas for sample {rep}...'.format(**locals()))
    start = configdir / ends.start
    end = configdir / ends.end
    rep_frame = gl.compute_gamma(start, end, controls, args.growth)
    rep_frame['rep'] = rep
    rep_frame = rep_frame.reset_index()
    rep_frames.append(rep_frame)
  collected = pd.concat(rep_frames).reset_index(drop=True)
  variants = collected.variant
  annoframe = gl.annotate_variants(variants, args.targetfile, args.locifile, args.genbank)
  geneorig = pd.DataFrame(annoframe[['gene', 'original']])
  geneorig = geneorig.dropna().drop_duplicates()
  annoframe = annoframe.drop('variant', axis='columns')
  annoframe = pd.concat([collected, annoframe], axis='columns')
  annoframe = annoframe.drop('original', axis='columns')
  origmap = pd.read_csv(args.origmap, sep='\t')
  annoframe = annoframe.merge(origmap, on='variant', how='left')
  annoframe = annoframe.dropna(subset=['original'])
  annoframe = annoframe.loc[double_mismatches(annoframe)]
  flatframe = flatgamma(annoframe, controls)
  flatframe = flatframe.drop(['locus_tag', 'gene'], axis='columns')
  flatframe = flatframe.reset_index()
  subvars = flatframe.apply(sub_variants, axis='columns')
  flatframe = pd.concat([flatframe, subvars], axis='columns')
  flatframe = flatframe.merge(geneorig, on='original', how='left')
  if args.justgene is not None:
    flatframe = flatframe.loc[flatframe.gene == args.justgene]
    assert len(flatframe) > 0
  aref = pd.DataFrame(flatframe[['original', 'vara']])
  aref.columns = ['original', 'variant']
  flatframe['a_pred'] = ml.predict_mismatch_scores(aref)
  bref = pd.DataFrame(flatframe[['original', 'varb']])
  bref.columns = ['original', 'variant']
  flatframe['b_pred'] = ml.predict_mismatch_scores(bref)
  flatframe['geometric'] = flatframe.a_pred * flatframe.b_pred
  flatframe['arithmetic'] = (flatframe.a_pred + flatframe.b_pred - 1).clip(0)
  flatframe['higher'] = flatframe[['a_pred', 'b_pred']].max(axis='columns')
  flatframe['lower'] = flatframe[['a_pred', 'b_pred']].min(axis='columns')
  flatframe['gap'] = flatframe.apply(bubble_gap, axis='columns')
  for compcol in ['geometric', 'arithmetic', 'higher', 'lower']:
    if args.justgene is not None:
      suffix = '.{args.justgene}.{compcol}.png'.format(**locals())
    else:
      suffix = '.{compcol}.png'.format(**locals())
    pngfile = pathlib.Path(args.pngfile).with_suffix(suffix)
    plot_gvk_doubles(flatframe, compcol, pngfile)

##############################################
if __name__ == "__main__":
  sys.exit(main())

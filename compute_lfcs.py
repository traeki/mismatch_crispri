#!/usr/bin/env python
# Author: John Hawkins (jsh) [really@gmail.com]

import argparse
import logging
import pathlib
import sys

import pandas as pd

import gamma_lib as gl
import model_lib as ml

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
      '--configdir', type=str, help='file: name of directory containing config.tsv',
      default=str(TESTDIR))
  parser.add_argument(
      '--lfcfile', type=str,
      help='file: file to which to write annotated log2-fold change measurements',
      default=None)
  parser.add_argument(
      '--growth', type=int,
      help='int: number of generations grown (in other words, g*t)',
      default=10)
  args = parser.parse_args()
  # TODO(jsh): Add check that either all or none of these are specified
  if args.lfcfile is None:
    args.lfcfile = pathlib.Path(args.configdir) / 'lfcs.tsv'
  return args


def flatlfc(stacked_replicates, controls):
  data = stacked_replicates
  data['y_pred'] = ml.predict_mismatch_scores(data)
  data.y_pred = data.y_pred.where(~data.variant.isin(controls), 0.0)
  anno = data.drop(['rep', 'lfc', 'start_mask'], axis='columns')
  anno = anno.drop_duplicates()
  lfc = data[['variant', 'rep', 'lfc']].set_index(['variant', 'rep'])
  lfc = lfc.unstack(level='rep').mean(axis='columns')
  data = pd.DataFrame(anno)
  data.set_index('variant', inplace=True)
  data['lfc'] = lfc
  return data


def count_mismatches(row):
  n = 0
  v = row.variant
  o = row.original
  for i in range(len(v)):
    if v[i] != o[i]:
      n += 1
  return n

def one_or_fewer_mismatches(row):
  return count_mismatches(row) <= 1

def max_one_mismatch_mask(frame):
  return frame.apply(one_or_fewer_mismatches, axis='columns')


def main():
  args = parse_args()
  controls = set(pd.read_csv(args.controls, header=None)[0])
  configdir = pathlib.Path(args.configdir)
  config = pd.read_csv(configdir / 'config.tsv', sep='\t')
  config = config.set_index('sample')
  rep_frames = list()
  for rep, ends in config.iterrows():
    logging.info('Computing lfcs for sample {rep}...'.format(**locals()))
    start = configdir / ends.start
    end = configdir / ends.end
    rep_frame = gl.compute_logdiff(start, end, controls)
    rep_frame['rep'] = rep
    rep_frame = rep_frame.reset_index()
    rep_frames.append(rep_frame)
  collected = pd.concat(rep_frames).reset_index(drop=True)
  variants = collected.variant
  annoframe = gl.annotate_variants(variants, args.targetfile, args.locifile, args.genbank)
  annoframe = annoframe.drop('variant', axis='columns')
  annoframe = pd.concat([collected, annoframe], axis='columns')
  annoframe.to_csv(args.lfcfile, index=False, sep='\t')
  annoframe = annoframe.loc[max_one_mismatch_mask(annoframe)]
  flatframe = flatlfc(annoframe, controls)
  flatfile = pathlib.Path(args.lfcfile).with_suffix('.mean.tsv')
  flatframe.to_csv(flatfile, sep='\t')

##############################################
if __name__ == "__main__":
  sys.exit(main())

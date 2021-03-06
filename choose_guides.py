#!/usr/bin/env python
# Author: John Hawkins (jsh) [really@gmail.com]

import argparse
import logging
import pathlib
import sys

import pandas as pd
import ipdb

import choice_lib as cl
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
      '--targetfile', type=str,
      help='file: ???.targets.all.tsv (as generated by traeki/sgrna_design)',
      default=str(TESTDIR / 'test.gb.targets.all.tsv'))
  parser.add_argument(
      '--locifile', type=str,
      help='file: list of locus_tag entries for which to design guides',
      default=str(TESTDIR / 'test.loci'))
  parser.add_argument(
      '--families', type=int,
      help='int: # of families for which to design variants (when available)',
      default=10)
  parser.add_argument(
      '--n', type=int,
      help='int: total # of variants to design',
      default=100)
  parser.add_argument(
      '--divide_evenly', action='store_true',
      help='Distribute guides within each family instead of each locus.')
  parser.add_argument(
      '--outfile', type=str,
      help='Distribute guides within each family instead of each locus.',
      default=str(TESTDIR / 'test.chosen.guides.tsv'))
  args = parser.parse_args()
  if args.divide_evenly and (args.n % args.families > 0):
    err = '--families:{args.families} not evenly divisible by --n:{args.n}'
    logging.warn(err.format(**locals()))
    guides_per_parent = args.n // args.families
    total = guides_per_parent * args.families
    template = '...reducing --n to {total} and preceeding anyway...'
    logging.warn(template.format(**locals()))
  return args

def main():
  args = parse_args()
  logging.info('Reading targets from {args.targetfile}...'.format(**locals()))
  logging.info('Building variants for {args.locifile}...'.format(**locals()))
  loci = set(pd.read_csv(args.locifile, sep='\t', header=None)[0])
  targetframe = pd.read_csv(args.targetfile, sep='\t')
  filtered = cl.filter_targets(targetframe, loci)
  pair_frame = cl.build_pairs(filtered, loci)
  pair_frame['y_pred'] = ml.predict_mismatch_scores(pair_frame)
  all_targets = pd.read_csv(args.targetfile, sep='\t')
  important = set(pd.read_csv(args.locifile, sep='\t', header=None)[0])
  chosen_loci = important
  # loop over locus tags and choose measure
  guides = dict()
  for locus in chosen_loci:
    template = 'Examining options for locus_tag: {locus}...'
    logging.info(template.format(**locals()))
    locus_preds = pair_frame.loc[pair_frame.locus_tag == locus]
    if len(locus_preds) == 0:
      logging.warn('...NO OPTIONS FOUND.')
      continue
    locus_targets = all_targets.loc[all_targets.locus_tag == locus]
    parents = cl.pick_n_parents(locus_preds, locus_targets, args.families)
    if args.divide_evenly:
      guides_per_parent = args.n // args.families
      guides[locus] = cl.choose_n_for_each(parents, locus_preds, guides_per_parent)
    else:
      parents = list(parents)
      candidates = locus_preds.loc[locus_preds.original.isin(parents)]
      guides[locus] = cl.choose_n_by_pred(candidates, args.n)
  allguides = set()
  for locus in guides:
    allguides.update(guides[locus])
  outframe = pair_frame.loc[pair_frame.variant.isin(allguides)]
  outframe.to_csv(args.outfile, sep='\t', index=False)

##############################################
if __name__ == "__main__":
  sys.exit(main())

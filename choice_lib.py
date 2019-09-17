#!/usr/bin/env python
# Author: John Hawkins (jsh) [really@gmail.com]

from collections import Counter
import logging
import random
import sys

import numpy as np
import pandas as pd

BASES = 'ACGT'
BIN_MIN = 0.1
BIN_MAX = 0.9
NBINS = 5

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s')
np.set_printoptions(precision=4, suppress=True)

def pred_bins():
  return np.linspace(BIN_MIN, BIN_MAX, NBINS-1)

def bin_preds(preds, bins):
  bins = np.digitize(preds, bins)
  return bins

def pick_n_parents(locus_preds, locus_targets, n):
  chosen = Counter()
  offset_map = locus_targets[['target', 'offset']]
  offset_map.columns = ['original', 'offset']
  offset_map.set_index('original', inplace=True)
  fullset = set(locus_preds.original.unique())
  past_picks = set()
  unused_picks = fullset - past_picks
  fallback_picks = set()
  # auto-include past picks but skip pre-filtered
  chosen.update(past_picks.intersection(fullset))
  while sum(chosen.values()) < n:
    used_offsets = offset_map.loc[offset_map.index.intersection(chosen.keys())]
    mask = offset_map.index.intersection(unused_picks)
    unused_offsets = offset_map.loc[mask]
    unused_offsets = unused_offsets.sort_values('offset')
    if unused_picks:
      candidate = unused_offsets.iloc[0]
      unused_picks.remove(candidate.name)
      unused_offsets = unused_offsets.iloc[1:]
      offset_diffs = used_offsets.offset - candidate.offset
      if (offset_diffs.abs() < 20).any():
        fallback_picks.add(candidate.name)
      else:
        chosen[candidate.name] += 1
    elif fallback_picks:
      pick = random.sample(fallback_picks, 1)[0]
      chosen[pick] += 1
      fallback_picks.remove(pick)
    else:
      ceiling = chosen.most_common(1)[0][-1] # most common count
      floor = chosen.most_common()[-1][-1] # least common count
      while True:
        # pick elements until we find one with sub-cap count (or we're level)
        boost = random.choice(chosen.most_common())
        if (floor == ceiling) or (boost[1] < ceiling):
          chosen[boost[0]] += 1
          break
  if chosen.most_common(1)[0][-1] > 4:
    # 5*10 + 13 > 60, so warn if we see 5+ instances
    template = 'Had to fall back to same family too often: {chosen}'
    logging.warn(template.format(**locals()))
  for ele in chosen.elements():
    yield ele

def choose_n_for_each(parents, preds, n):
  chosen = set()
  for parent in parents:
    needed = n
    family_preds = preds.loc[preds.original == parent]
    used_mask = family_preds.variant.isin(chosen)
    family_remaining = family_preds.loc[~used_mask]
    picks = choose_n_by_pred(family_remaining, needed)
    chosen.update(picks)
  return chosen

def choose_n_by_pred(preds, n):
  return choose_n_by_bin(preds, 'y_pred', n)

def choose_n_by_bin(data, binnable, n):
  if data.empty:
    return set()
  loci = set(data.locus_tag.unique())
  locus = list(loci)[0]
  if len(loci) != 1:
    template = 'choose_n_by_bin called with multiple loci: {loci}'
    logging.fatal(template.format(**locals()))
    sys.exit(2)
  usable = data.dropna()
  usable = usable.copy()
  if usable.shape[0] < n:
    n_found = usable.shape[0]
    template = 'Only found {n_found}/{n} binnable guides for locus {locus}'
    logging.warn(template.format(**locals()))
    if data.shape[0] < n:
      template = 'Fewer than {n} guides exist for locus {locus}'
      logging.warning(template.format(**locals()))
    return set(usable.variant)
  # ascribe bins
  bins = pred_bins()
  usable['bin'] = bin_preds(usable[binnable], bins)
  # choose guides for each bin (skipping 0)
  chosen = set()
  bins = list(range(NBINS))
  per_bin = (n-1) // (len(bins)-1)
  for b in bins:
    k = per_bin
    if b == bins[-1]:
      k -= 1
    bin_items = usable.loc[usable.bin == b]
    if bin_items.shape[0] <= k:
      # if there aren't more than k items in this bin, grab what we have
      chosen.update(bin_items.variant)
    else:
      # if there are more than k, pick at random
      poss = set(bin_items.variant) - chosen
      chosen.update(random.sample(poss, k))
  # How many more do we need?
  z = (n - len(chosen))
  # Grab up to z preferring non-max efficacy
  leftover = set(usable.variant) - chosen
  toosick = set(usable.loc[usable.bin == bins[-1]].variant)
  okset = leftover - toosick
  if len(okset) >= z:
    chosen.update(random.sample(okset, z))
  else:
    chosen.update(okset)
    dregs = toosick - chosen
    chosen.update(random.sample(dregs, (z - len(okset))))
  assert len(chosen) == n
  return chosen


def all_single_variants(parents):
  pairs = list()
  for parent in parents:
    assert len(parent) == 20
    children = set()
    for i in range(len(parent)):
      for letter in BASES:
        children.add(parent[:i] + letter + parent[i+1:])
    children.remove(parent)
    for child in children:
      pairs.append((parent, child))
  pairrows = list()
  for original, variant in pairs:
    pairrows.append({'original':original, 'variant':variant})
  return pd.DataFrame(pairrows)

def filter_targets(parentframe, loci):
  antisense_rows = parentframe.loc[parentframe.transdir=='anti']
  important_rows = antisense_rows.loc[antisense_rows.locus_tag.isin(loci)]
  dupmask = ~important_rows.target.duplicated(keep=False)
  unique_rows = important_rows.loc[dupmask]
  return unique_rows

def build_pairs(targetframe, loci):
  parents = targetframe.target
  origvars = all_single_variants(parents)
  locmap = targetframe[['target', 'locus_tag']].set_index('target').locus_tag
  pammap = targetframe[['target', 'pam']].set_index('target').pam
  origvars['locus_tag'] = origvars.original.map(locmap)
  origvars['pam'] = origvars.original.map(pammap)
  return origvars

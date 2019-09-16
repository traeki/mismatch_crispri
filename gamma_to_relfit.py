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
      '--gammafile', type=str,
      help='file: file with gamma scores',
      required=True)
  parser.add_argument(
      '--relfitfile', type=str,
      help='file: destination for relfit scores',
      required=True)
  args = parser.parse_args()
  return args


def main():
  args = parse_args()
  gammaframe = pd.read_csv(args.gammafile, sep='\t')
  relfitframe = gammaframe.copy(deep=True)
  relfitframe.gamma = relfitframe.gamma + 1
  relfitframe.rename({'gamma':'relfit'}, axis='columns', inplace=True)
  relfitframe.to_csv(args.relfitfile, sep='\t', index=False)

##############################################
if __name__ == "__main__":
  sys.exit(main())

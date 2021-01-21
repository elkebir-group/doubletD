#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 5 2020

@author: Palash Sashittal
"""

import pandas as pd
import scrublet as scr
import numpy as np
import sys
import argparse
import itertools
import math

def main(args):

    df = pd.read_csv(args.inputTernary, sep='\t')
    genotype_mat = df[df.columns[1:]].values.transpose() 
    app_genotype_mat = np.hstack((np.array(range(genotype_mat.shape[0]))[:, np.newaxis], genotype_mat)) 
    np.savetxt(args.outputfile, app_genotype_mat, fmt='%d')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputTernary", type=str, help="csv file with a table of alternate read counts for each position in each cell")
    parser.add_argument("-o", "--outputfile", type=str, help="output file name")

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)

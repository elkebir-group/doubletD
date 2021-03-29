#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 2020

@author: Palash Sashittal
"""

import pandas as pd
import sys
import argparse
import itertools
import math
import numpy as np

class doubletFinder():

    def __init__(self, df_total, df_alt, delta, beta, mu_hetero, mu_hom, alpha_fp, alpha_fn, missing=False, verbose=True, binom=False, precision=None):

        self.df_total = df_total
        self.df_alt = df_alt
        self.delta = delta
        self.beta = beta
        self.mu_hetero = mu_hetero
        self.mu_hom = mu_hom
        self.alpha_fp = alpha_fp
        self.alpha_fn = alpha_fn
        self.precision = precision
        self.missing = missing
       
        if binom:
            self.prv_y = self.prv_y_b
        else:
            self.prv_y = self.prv_y_bb

        self.cells = list(df_total['cell_id'].values)
        self.muts = list(df_total.columns[1:])
        self.df_total = self.df_total.set_index(['cell_id'])
        self.df_alt = self.df_alt.set_index(['cell_id'])

        if (self.df_total.values - self.df_alt.values).min() < 0:
            raise Exception('total reads must be greater than or equal to alternate reads!')

        self.Sigma = ((0, 1/2, 1), (0, 1/4, 1/2, 3/4, 1))
        self.Theta = ((0, 1/2, 1), (0, 1/4, 1/3, 1/2, 2/3, 3/4, 1))

        self.px_z = {x: 0 for z in [0,1] for x in itertools.product(self.muts, self.Sigma[z], [z])}

        if self.alpha_fn == None or self.alpha_fp == None or self.precision == None:
            # get vaf numpy array
            vaf_values = (self.df_alt / self.df_total).values
            vaf_values = vaf_values[~np.isnan(vaf_values)]

            # filter vafs less than 0.15 and get alpha_fp
            if self.alpha_fp == None or self.precision == None:
                mean1, prec1 = getBetaMOM(vaf_values[vaf_values <= 0.15])
                self.alpha_fp = mean1
            
            if self.alpha_fn == None or self.precision == None:
                mean2, prec2 = getBetaMOM(vaf_values[vaf_values >= 0.85])
                self.alpha_fn = 1 - mean2

            if self.precision == None:
                mean3, prec3 = getBetaMOM(vaf_values[(vaf_values > 0.15) & (vaf_values < 0.85)])
                self.precision = np.median([prec1, prec2, prec3])

        print(f"alpha_fn = {self.alpha_fn}")
        print(f"alpha_fp = {self.alpha_fp}")
        print(f"precision = {self.precision}")

        if self.mu_hom == None or self.mu_hetero == None:
            estimate = True
        else:
            estimate = False

        for m in self.muts:
            if estimate:
            
                reads = self.df_alt[m]
                total = self.df_total[m]
                vaf = pd.DataFrame(reads/total)
            
                loh_cells = vaf.loc[vaf[m] > 0.85]
                wt_cells = vaf.loc[vaf[m] < 0.15]
                het_cells = vaf.loc[(vaf[m] >= 0.15 )& (vaf[m] <= 0.85)]

                
                # #count the total number of non-na VAF cells for variant m
                non_na_cells =  loh_cells.shape[0] + wt_cells.shape[0] + het_cells.shape[0]
                est_wt_rate = wt_cells.shape[0]/non_na_cells
                est_loh_rate = loh_cells.shape[0]/non_na_cells
                est_het_rate = het_cells.shape[0]/non_na_cells
       
                self.px_z[m, 0,0] = est_wt_rate
                self.px_z[m, 1/2, 0] =  est_het_rate
                self.px_z[m, 1, 0] = est_loh_rate


            else:
                self.px_z[m, 0,0] = 1 - self.mu_hom - self.mu_hetero
                self.px_z[m, 1/2,0] =  self.mu_hetero
                self.px_z[m, 1, 0] = self.mu_hom
            

        norm_const = {}
        for m in self.muts:
            for a,b in itertools.product(self.Sigma[0], repeat = 2):
                c = (a + b)/2
                if c in self.Sigma[1]:
                    self.px_z[m,c,1] += self.px_z[m,a,0] * self.px_z[m,b,0]
            norm_const[m] = sum([self.px_z[m,a,0] * self.px_z[m,b,0] for a,b in itertools.product(self.Sigma[0], repeat = 2)])

        for m in self.muts:
            for c in self.Sigma[1]:
                self.px_z[m,c,1] /= norm_const[m]
        

        self.py_xz = {x: 0 for z in [0,1] for x in itertools.product(self.Theta[z], self.Sigma[z], [z])}

        self.py_xz[0,   0,   0]= 1 - self.beta**2
        self.py_xz[0,   1/2, 0]= self.beta * (1 - self.beta)
        self.py_xz[1/2, 1/2, 0]= (1-self.beta)**2
        self.py_xz[1,   1/2, 0]= self.beta * (1 - self.beta)
        self.py_xz[1,   1,   0]= 1 - self.beta**2
	
        self.py_xz[0,   0,   1]= 1 - self.beta**4
        self.py_xz[0, 1/4,   1]= self.beta * (1 - self.beta)**3 + 3 * self.beta**2 * (1 - self.beta)**2 + 3 * self.beta**3 * (1 - self.beta)
        self.py_xz[1/4, 1/4, 1]= (1 - self.beta)**4 
        self.py_xz[1/3, 1/4, 1]= 3*self.beta*(1 - self.beta)**3
        self.py_xz[1/2, 1/4, 1]= 3*self.beta**2 * (1 - self.beta)**2
        self.py_xz[1,   1/4, 1]= self.beta**3 * (1 - self.beta)
        self.py_xz[1/3, 1/2, 1]= 2 * self.beta * (1 - self.beta)**3
        self.py_xz[1/2, 1/2, 1]= (1-self.beta)**4 + 4 * self.beta**2 * (1 - self.beta)**2
        self.py_xz[2/3, 1/2, 1]= 2 * self.beta * (1 - self.beta)**3
        self.py_xz[1,   1/2, 1]= self.beta**2 * (1 - self.beta)**2 + 2 * self.beta**3 * (1 - self.beta)
        self.py_xz[0,   3/4, 1]= self.py_xz[1,   1/4, 1]
        self.py_xz[1/2, 3/4, 1]= self.py_xz[1/2, 1/4, 1]
        self.py_xz[2/3, 3/4, 1]= self.py_xz[1/3, 1/4, 1]
        self.py_xz[3/4, 3/4, 1]= self.py_xz[1/4, 1/4, 1]
        self.py_xz[1,   3/4, 1]= self.py_xz[0,   1/4, 1]
        self.py_xz[1, 1, 1] = 1 - self.beta**4



        self.doublet_result = None

        if self.delta == 0:
            self.threshold = sys.float_info.max
        elif self.delta == 1:
            self.threshold = -sys.float_info.max
        else:
            self.threshold = math.log((1 - self.delta) / self.delta)

    def solve(self):

        self.doublet_result = {}
        self.logprobs = {}
        for cell in self.cells:
            self.logprobs[cell, 0] = self.prv_z(cell, 0)
            self.logprobs[cell, 1] = self.prv_z(cell, 1)

            if self.logprobs[cell, 1] - self.logprobs[cell, 0] > self.threshold:
                self.doublet_result[cell] =  'doublet'
            else:
                self.doublet_result[cell] = 'singlet'

    def prv_z(self, cell, z):

        log_prob_sum = 0
        for mut in self.muts:
            v = self.df_alt.loc[cell, mut]
            r = self.df_total.loc[cell, mut] - v
            if self.missing and r + v == 0:
                if z == 0 and self.beta > 0:
                    log_prob_sum += 2 * math.log(self.beta)
                if z == 1 and self.beta > 0:
                    log_prob_sum += 4 * math.log(self.beta)
                continue
            prob_sum = 0
            for x in self.Sigma[z]:
                prob_sum_x = 0
                for y in self.Theta[z]:
                    prob_sum_x += self.prv_y(r, v, y) * self.py_xz[y, x, z]
                prob_sum += self.px_z[mut, x, z] * prob_sum_x

            log_prob_sum += math.log(prob_sum)
        return log_prob_sum

    def prv_y_b(self, r, v, y):
        yprime = self.alpha_fp + (1 - self.alpha_fp - self.alpha_fn) * y
        return nCr(r+v, v) * (yprime ** v) * ((1-yprime) ** r)

    def prv_y_bb(self, r, v, y):
        yprime = self.alpha_fp + (1 - self.alpha_fp - self.alpha_fn) * y
        if yprime == 0:
            yprime = 0.001
        if yprime == 1:
            yprime = 0.999
        alpha = self.precision*yprime 
        beta = self.precision - alpha
        n = r + v 
        #print(f"n:{n} r:{r} v:{v} p:{y} alpha:{alpha} beta:{beta}")
        num = math.lgamma(n+1) + math.lgamma(v+alpha) + math.lgamma(n-v+beta) + math.lgamma(alpha+beta)
        den = math.lgamma(v+1) + math.lgamma(n-v+1) + math.lgamma(n+ alpha + beta) + math.lgamma(alpha) + math.lgamma(beta)
        prob = math.exp(num- den)

        return prob

    def writeSolution(self, outputFile):

        with open(outputFile, 'w') as output:
            output.write("cell_id\tprob_z0\tprob_z1\tprediction\n")
            for cell in self.cells:
                output.write(f"{cell}\t{self.logprobs[cell, 0]}\t{self.logprobs[cell,1]}\t{self.doublet_result[cell]}\n")

    def likelihood(self):
        likelihood = len(self.cells) * math.log(1 - self.delta) 
        for cell in self.cells:
            if self.doublet_result[cell] == 'doublet':
                likelihood += self.logprobs[cell,1] - self.logprobs[cell,0] - self.threshold
        return likelihood

def getBetaMOM(x):

    m_x = np.mean(x)
    s_x = np.std(x)
    x_alpha = m_x*((m_x*(1 - m_x)/s_x**2) - 1)
    x_beta = (1 - m_x)*((m_x*(1 - m_x)/s_x**2) - 1) 
    
    return x_alpha/(x_alpha + x_beta), x_alpha + x_beta

def nCr(n,r):
    f = math.factorial
    #print(n, r)
    return f(n) // f(r) // f(n-r)

def main(args):

    df_total = pd.read_csv(args.inputTotal)
    df_alt = pd.read_csv(args.inputAlternate)

    if len(df_total) != len(df_alt):
        raise Exception("number of cells in the two input files do not match!")

    if len(df_total.columns) != len(df_alt.columns):
        raise Exception("number of cells in the two input files do not match!")

    ncells = len(df_total)
    npos = len(df_total.columns) - 1

    if args.verbose:
        print(f"number of cells is {ncells}")
        print(f"number of mutation positions is {npos}")

    solver = doubletFinder(df_total, df_alt, delta = args.delta, beta = args.beta, missing = args.missing, 
                           mu_hetero = args.mu_hetero, mu_hom = args.mu_hom, alpha_fp = args.alpha_fp, alpha_fn = args.alpha_fn,
                           verbose = args.verbose, binom = args.binom, precision = args.prec)

    solver.solve()

    solver.writeSolution(args.outputfile)

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputTotal", type=str, help="csv file with a table of total read counts for each position in each cell")
    parser.add_argument("--inputAlternate", type=str, help="csv file with a table of alternate read counts for each position in each cell")
    parser.add_argument("--delta", type=float, default=0.1, help="expected doublet rate [0.1]")
    parser.add_argument("--beta", type=float, default=0.05, help="allelic dropout (ADO) rate [0.05]")
    parser.add_argument("--mu_hetero", type=float, help="heterozygous mutation rate [None]")
    parser.add_argument("--mu_hom", type=float, help="homozygous mutation rate [None]")
    parser.add_argument("--alpha_fp", type=float, help="copy false positive error rate [None]")
    parser.add_argument("--alpha_fn", type=float, help="copy false negative error rate [None]")
    parser.add_argument("-o", "--outputfile", type=str, help="output file name")
    parser.add_argument("--noverbose", dest="verbose", help="do not output statements from internal solvers [default is false]", action='store_false')
    parser.add_argument("--binomial", dest="binom", help="use binomial distribution for read count model [default is false]", action='store_true')
    parser.add_argument("--prec", type=float, help="precision for beta-binomial distribution [None]")
    parser.add_argument("--missing", dest="missing", help="use missing data in the model? [No]", action = 'store_true')
    parser.set_defaults(missing=False)
    parser.set_defaults(binom=False)
    parser.set_defaults(verbose=True)

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)

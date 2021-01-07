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
import scipy.special as ss

class doubletFinder():

    def __init__(self, df_total, df_alt, delta, beta, mu_hetero, mu_homo, alpha_fp, alpha_fn, asym=False, all_imb=False, missing=False, verbose=True, cellcoal=False, binom=False, precision=10000, estimate=False):

        self.df_total = df_total
        self.df_alt = df_alt
        self.delta = delta
        self.beta = beta
        self.mu_hetero = mu_hetero
        self.mu_homo = mu_homo
        self.alpha_fp = alpha_fp
        self.alpha_fn = alpha_fn
        self.precision = precision
        self.asym = asym
        self.missing = missing
       
        if binom:
            self.prv_y = self.prv_y_b
        else:
            if all_imb:
                self.prv_y = self.prv_y_all_imb
            else:
                self.prv_y = self.prv_y_bb

        self.cells = list(df_total['cell_id'].values)
        self.muts = list(df_total.columns[1:])
        self.df_total = self.df_total.set_index(['cell_id'])
        self.df_alt = self.df_alt.set_index(['cell_id'])

        self.Sigma = ((0, 1/2, 1), (0, 1/4, 1/2, 3/4, 1))
        self.Theta = ((0, 1/2, 1), (0, 1/4, 1/3, 1/2, 2/3, 3/4, 1))

        self.px_z = {x: 0 for z in [0,1] for x in itertools.product(self.muts, self.Sigma[z], [z])}
        if estimate:
            
            for m in self.muts:
                if estimate:
                
                    reads = self.df_alt[m]
                    total = self.df_total[m]
                    vaf = pd.DataFrame(reads/total)
                
                    loh_cells = (1-self.beta)*vaf.loc[vaf[m] > 0.85]
                    wt_cells = (1-self.beta)*vaf.loc[vaf[m] < 0.15]

                    #count the total number of non-na VAF cells for variant m
                    non_na_cells =  len(self.cells) -vaf[m].isna().sum()
                
                    est_wt_rate = wt_cells.shape[0]/non_na_cells
                    est_loh_rate = loh_cells.shape[0]/non_na_cells
                    
                    #print(f"variant:{m} wt:{est_wt_rate} het:{1-est_wt_rate - est_loh_rate} loh:{est_loh_rate}")
                    self.px_z[m, 0,0] = est_wt_rate
                    self.px_z[m, 1/2,0] =  1 -est_wt_rate - est_loh_rate
                    self.px_z[m, 1, 0] = est_loh_rate


                else:
                    self.px_z[m, 0,0] = 1 - self.mu_homo - self.mu_hetero
                    self.px_z[m, 1/2,0] =  self.mu_hetero
                    self.px_z[m, 1, 0] = self.mu_homo
            

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

        if not cellcoal:
            self.py_xz[0,   0,   0]= 1
            self.py_xz[0,   1/2, 0]= self.beta/2
            self.py_xz[1/2, 1/2, 0]= 1-self.beta
            self.py_xz[1,   1/2, 0]= self.beta/2
            self.py_xz[1,   1,   0]= 1

            self.py_xz[0,   0,   1]= 1
            self.py_xz[0, 1/4,   1]= self.beta/4
            self.py_xz[1/4, 1/4, 1]= 1-self.beta
            self.py_xz[1/3, 1/4, 1]= 3*self.beta/4
            self.py_xz[1/3, 1/2, 1]= self.beta/2
            self.py_xz[1/2, 1/2, 1]= 1-self.beta
            self.py_xz[2/3, 1/2, 1]= self.beta/2
            self.py_xz[2/3, 3/4, 1]= 3*self.beta/4
            self.py_xz[3/4, 3/4, 1]= 1-self.beta
            self.py_xz[1,   3/4, 1]= self.beta/4
            self.py_xz[1, 1, 1] = 1
        else:
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


        #print(self.py_xz)

        # self.px_z = {x: 0 for z in [0,1] for x in itertools.product(self.Sigma[z], [z])}

        # self.px_z[0,0] = 1 - self.mu_homo - self.mu_hetero
        # self.px_z[1/2,0] = self.mu_hetero
        # self.px_z[1, 0] = self.mu_homo

        # for a,b in itertools.product(self.Sigma[0], repeat = 2):
        #     c = (a + b)/2
        #     if c in self.Sigma[1]:
        #         self.px_z[c,1] += self.px_z[a,0] * self.px_z[b,0]
        # norm_const = sum([self.px_z[a,0] * self.px_z[b,0] for a,b in itertools.product(self.Sigma[0], repeat = 2)])

        # for c in self.Sigma[1]:
        #     self.px_z[c,1] /= norm_const

        #print(self.px_z)

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
            #print(f"{cell}\t{mut}\t{r}\t{v}\t{prob_sum}")
            log_prob_sum += math.log(prob_sum)
        return log_prob_sum

    def prv_y_b(self, r, v, y):
#        if self.asym:
#            yprime = y - 4 * self.alpha * y + 3 * self.alpha
#        else:
#            yprime = y - 2 * self.alpha * y + self.alpha
        yprime = self.alpha_fp + (1 - alpha_fp - alpha_fn) * y
        return nCr(r+v, v) * (yprime ** v) * ((1-yprime) ** r)

    def prv_y_bb(self, r, v, y):
#        if self.asym:
#            yprime = y - 4 * self.alpha * y + 3 * self.alpha
#        else:
#            yprime = y - 2 * self.alpha * y + self.alpha
        yprime = self.alpha_fp + (1 - alpha_fp - alpha_fn) * y
        if yprime == 0:
            yprime = 0.001
        alpha = self.precision*yprime 
        beta = self.precision - alpha
        n = r + v 
        #print(f"n:{n} p:{y} alpha:{alpha} beta:{beta}")
        num = math.lgamma(n+1) + math.lgamma(v+alpha) + math.lgamma(n-v+beta) + math.lgamma(alpha+beta)
        den = math.lgamma(v+1) + math.lgamma(n-v+1) + math.lgamma(n+ alpha + beta) + math.lgamma(alpha) + math.lgamma(beta)
        prob = math.exp(num- den)

        return prob

    def prv_y_bb_rev(self, r, v, y):
        if self.asym:
            a = 3 * self.alpha
            b = 1 - 4 * self.alpha
        else:
            a = self.alpha
            b = 1 - 2 * self.alpha

        alpha = self.precision * y
        beta = self.precision - alpha

        k = np.arange(v)
        l = np.arange(r)

        k_terms = ss.gamma(alpha + k) * ss.comb(v, k) * a**(v - k)
        l_terms = ss.gamma(beta + l) * ss.comb(r, l) * (1 - a - b)**(r - l)
        k = k.reshape(-1, 1)
        l = l.reshape(1, -1)
        kl_terms= b**(k + l) / ss.gamma(alpha + beta + k + l)

        return np.einsum("i,j,ij", k_terms, l_terms, kl_terms) * ss.comb(n, v) / ss.beta(alpha, beta)


    def prv_y_all_imb(self, r, v, y):
        if y in [0, 1]:
            return self.prv_y_b(r, v, y)
        else:
            return self.prv_y_bb(r, v, y)

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

    solver = doubletFinder(df_total, df_alt, delta = args.delta, beta = args.beta, all_imb = args.all_imb, missing = args.missing, 
                           mu_hetero = args.mu_hetero, mu_homo = args.mu_homo, alpha_fp = args.alpha_fp, alpha_fn = args.alpha_fn, asym = args.asym,
                           verbose = args.verbose, cellcoal = args.cellcoal, binom = args.binom, precision = args.prec,
                           estimate = args.estimate)

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
    parser.add_argument("--delta", type=float, default=0.1, help="doublet rate [0.1]")
    parser.add_argument("--beta", type=float, default=0.05, help="Allelic dropout (ADO) rate [0.05]")
    parser.add_argument("--mu_hetero", type=float, default=0.5, help="heterozygous mutation rate [0.5]")
    parser.add_argument("--mu_homo", type=float, default=0, help="homozygous mutation rate [0]")
    parser.add_argument("--alpha_fp", type=float, default = 0, help="copy false positive error rate [0]")
    parser.add_argument("--alpha_fn", type=float, default = 0, help="copy flase negative error rate [0]")
    parser.add_argument("-o", "--outputfile", type=str, help="output file name")
    parser.add_argument("--noverbose", dest="verbose", help="do not output statements from internal solvers [default is false]", action='store_false')
    parser.add_argument("--cellcoal", dest="cellcoal", help="use cellcoal doublet model [default is false]", action='store_true')
    parser.add_argument("--binomial", dest="binom", help="use cellcoal doublet model [default is false]", action='store_true')
    parser.add_argument("--prec", type=float, default = 10000, help="Precision for Beta Distribution [10000]")
    parser.add_argument("--asym", dest="asym", help="use asymmetric sequencing error model? [No]", action = 'store_true')
    parser.add_argument("--allelic-imbalance", dest="all_imb", help="allelic imbalance model? [No]", action = 'store_true')
    parser.add_argument("--missing", dest="missing", help="use missing data in the model? [No]", action = 'store_true')
    parser.add_argument("--estimate", dest="estimate", help="Estimate mutation rates from input data", action = 'store_true')
    parser.set_defaults(missing=False)
    parser.set_defaults(all_imb=False)
    parser.set_defaults(asym=False)
    parser.set_defaults(binom=False)
    parser.set_defaults(cellcoal=False)
    parser.set_defaults(estimate=False)
    parser.set_defaults(verbose=True)

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)

#!/usr/bin/env python

"""
Functions for creating a SeqLib Class object.
"""

import numpy as np
import pandas as pd


class Seqlib:
    def __init__(self, ninds, nsites):
        self.ninds = ninds
        self.nsites = nsites
        self.arr = self.simulate()

    def mutate(self, base):
        diff = set("ACTG") - set(base)# provides a set of 4 bases and subtracts an individual base from the original set 
        return np.random.choice(list(diff))# returns one of the three remaining bases after the base is subtracted
        
    def simulate(self):
        oseq = np.random.choice(list("ACGT"), size=self.nsites) # calls orginal sequence array and creates a sequence of bases
        self.arr = np.array([oseq for i in range(self.ninds)]) # creates new array of sequences for the range of # of individuals 
        muts = np.random.binomial(1, 0.1, (self.ninds, self.nsites)) # iterates over columns and creates mutation with a 10% mutation probability 
    
        for col in range(self.nsites):
            newbase = mutate(self.arr[0, col]) #Use mutate function to create mutations in first row of each column
            mask = muts[:, col].astype(bool)#goes over each column and grabs only the sites that are mutated - returns boolean type
            self.arr[:, col][mask] = newbase# grab the column, apply mask to .....
        missing = np.random.binomial(1, 0.1, (self.ninds, self.nsites)) # create missing values using a binomial distribution 
        self.arr[missing.astype(bool)] = "N" # return missing values with "N"
        return self.arr # return array

    def filter_missing(self, maxfreq):
        freqmissing = np.sum(self.arr == "N", axis=0) / self.arr.shape[0]#sum every first row element that equals "N", divide by the number of rows 
        return self.arr[:, freqmissing <= maxfreq] # return all rows in array where the proportion of freqmissing equals the value given for maxfreq

    def filter_maf(self, minfreq):
        freqs = np.sum(self.arr != self.arr[0], axis=0) / self.arr.shape[0]# sum every row element that doesn't match the first element in the column, divide by the number of rows 
        maf = freqs.copy() # make copy to not change original frequency array
        maf[maf > 0.5] = 1 - maf[maf > 0.5] # subselect sites with major freq (>0.5) and modify to be 1-value
        return self.arr[:, maf > minfreq] # return array rows where major allele frequency is greater than minimum allele frequency
    
    def descriptive(self, minmaf, maxmissing):
        job1 = self.filter_missing(maxfreq = maxmissing)
        freqs = np.sum(job1 != job1[0], axis=0) / job1.shape[0]  
        maf = freqs.copy() 
        maf[maf > 0.5] = 1 - maf[maf > 0.5] 
        return job1[:, maf > minmaf] 
    
    def calculcate_statistics(self):
        nd = np.var(self.arr == self.arr[0], axis=0).mean() # get mean value for all elements present
        mf = np.mean(np.sum(self.arr != self.arr[0], axis=0) / self.arr.shape[0]) # get mean value for elements that have mutated (=1)
        inv = np.any(self.arr != self.arr[0], axis=0).sum() # get sum value for all elements that have not mutated (=0)
        var = self.arr.shape[1] - inv # subtract invariant sites from whole array to get the variable sites
        return pd.Series(
        {"mean nucleotide diversity": nd,
         "mean minor allele frequency": mf,
         "invariant sites": inv,
         "variable sites": var,
        })
    
    
    
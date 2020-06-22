# -*- coding: utf-8 -*-
"""
    Example data set for the sgpv-package
    This is the same dataset that is used in the R-library of the same name
Created on Tue Jun 16 16:20:27 2020

@author: Kristjan
 Test Statistics from Gloub (1999) Leukemia data set

 Data are from 7218 gene specific t-tests for a difference in mean expression 
 (on the log scale; AML versus ALL) in the Gloub data set (1999).
 Data are from 72 patients using a pooled t-test (df=70). 
 Included in the dataframe are the following:
     t-statistic (t.stat), p-value (p.value), CI lower limit (ci.lo),
     CI upper limit (ci.hi), estimate (estimate), standard error (se).


 @usage data(leukstats)

 @format An object of class data.frame. Includes the following: 
     t-statistic (t.stat), p-value (p.value), 
     CI lower limit (ci.lo), CI upper limit (ci.hi), 
     estimate (estimate), standard error (se).

 @references Gloub (1999) and used in Blume et. al. (2018) PlosONE.

 Blume JD, D’Agostino McGowan L, Dupont WD, Greevy RA Jr. (2018). Second-generation p-values: Improved rigor, reproducibility, & transparency in statistical analyses. \emph{PLoS ONE 13(3): e0188299. https://doi.org/10.1371/journal.pone.0188299

 @source https://github.com/ramhiser/datamicroarray/wiki/Golub-(1999)
"""

import pandas as pd
df = pd.read_csv('../data/leukstats.csv', index_col=0)
df.sort_values('p.value')

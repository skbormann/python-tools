# -*- coding: utf-8 -*-
"""
Example data set for the sgpv-module.

This is the same dataset that is used in the R-library of the same name.
Test Statistics from Gloub (1999) Leukemia data set.

Data are from 7218 gene specific t-tests for a difference in mean expression
(on the log scale; AML versus ALL) in the Gloub data set (1999).
Data are from 72 patients using a pooled t-test (df=70).
Included in the dataframe are the following:
t-statistic (t.stat), p-value (p.value), CI lower limit (ci.lo),
CI upper limit (ci.hi), estimate (estimate), standard error (se).


References:
 Gloub (1999) and used in Blume et. al. (2018) PlosONE.

 Blume JD, D’Agostino McGowan L, Dupont WD, Greevy RA Jr. (2018).
 Second-generation p-values: Improved rigor, reproducibility, & transparency
 in statistical analyses. PLoS ONE 13(3): e0188299. https://doi.org/10.1371/journal.pone.0188299

 Source: https://github.com/ramhiser/datamicroarray/wiki/Golub-(1999)
"""
# Copied code from statsmodels.datasets.utils to add allow import of dataset from
# module directory

def load():
    """
    Load the example

    Returns
    -------
    df : pandas_series
        A dataframe containing the dataset.

    """
    df = load_csv(__file__, 'leukstats.csv', index_col=0, convert_float=True)
    return df


# Directly copied from statsmodels.datasets.utils and added an index_col argument
def load_csv(base_file, csv_name, index_col, sep=',', convert_float=False):
    from os.path import dirname, abspath, join
    from pandas import read_csv
    """Standard simple csv loader"""
    filepath = dirname(abspath(base_file))
    filename = join(filepath, csv_name)
    engine = 'python' if sep != ',' else 'c'
    float_precision = {}
    if engine == 'c':
        float_precision = {'float_precision': 'high'}
    data = read_csv(filename, sep=sep, engine=engine, index_col=index_col,
                    **float_precision)
    if convert_float:
        data = data.astype(float)
    return data

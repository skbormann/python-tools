# sgpv module

This module allows to calculate Second Generation P-Values developed by Blume et.al.(2018,2019) 
and their associated diagnostics in Python.
This package is a translation of the original [sgpv R-library](https://github.com/weltybiostat/sgpv) into Python.
The same library has already been translated into [Stata](https://github.com/skbormann/stata-tools/sgpv) by the author of this Python translation.

This module contains the following functions:

            value    - calculate the SGPVs
            power    - power functions for the SGPVs
            risk     - false confirmation/discovery risks for the SGPVs
            plot     - plot the SGPVs
            data     - load the example dataset into memory

The module comes with an example dataset (leukstats.csv) to showcase the plotting function.
See the documentation in the file data.py for more information about this dataset.

# Dependencies

This module depends on: 
 pandas>=1.0.4, matplotlib>=3.2.1, numpy>=1.18.0, scipy>=1.3.2
            
These dependencies document only under which version I tested my functions.
Older version might work as well. 

# Installation
Binaries and source distributions are available from PyPi
[https://pypi.org/projects/sgpv](https://pypi.org/projects/sgpv)


The same installation files are also located in the folder [dist](https://github.com/skbormann/python-tools/sgpv/dist).
Just download the tarball and unzip it. Then run
```python
python setup.py install
```

# Examples
Below are some examples taken from the documentation of each function:

## Calculate second generation p-values (sgpv.value):
```python
>>> import numpy as np
>>> from sgpv import sgpv
>>> lb = (np.log(1.05), np.log(1.3), np.log(0.97))
>>> ub = (np.log(1.8), np.log(1.8), np.log(1.02))
>>> sgpv.value(est_lo = lb, est_hi = ub,
             null_lo = np.log(1/1.1), null_hi = np.log(1.1))
    sgpv(pdelta=array([0.1220227, 0.        , 1.        ]),
     deltagap=array([None, 1.7527413, None], dtype=object))
```

## Power function (sgpv.power): 
```python
>>> from sgpv import sgpv       
>>> sgpv.power(true=2, null_lo=-1, null_hi=1, std_err = 1,
...        interval_type='confidence', interval_level=0.05)
poweralt = 0.168537 powerinc = 0.831463 powernull =  0
type I error summaries:
at 0 = 0.0030768 min = 0.0030768 max = 0.0250375 mean = 0.0094374
>>> sgpv.power(true=0, null_lo=-1, null_hi=1, std_err = 1,
...         interval_type='confidence', interval_level=0.05)
poweralt = 0.0030768 powerinc = 0.9969232 powernull =  0
type I error summaries:
at 0 = 0.0030768 min = 0.0030768 max = 0.0250375 mean = 0.0094374 
```

## False discory risk (sgpv.risk):   
```python           
>>> from sgpv import sgpv
>>> import numpy as np
>>> from scipy.stats import norm
>>> sgpv.risk(sgpval = 0, null_lo = np.log(1/1.1), null_hi = np.log(1.1),
           std_err = 0.8, null_weights = 'Uniform',
           null_space = (np.log(1/1.1), np.log(1.1)), alt_weights = 'Uniform',
           alt_space = (2 + 1*norm.ppf(1-0.05/2)*0.8, 2 - 1*norm.ppf(1-0.05/2)*0.8),
           interval_type = 'confidence', interval_level = 0.05);
The false discovery risk (fdr) is: 0.0594986
```

## Plotting of SGPVs with example dataset (sgpv.plot):
```python
>>> from sgpv import sgpv
>>> from sgpv import data
>>> import matplotlib.pyplot as plt
>>> df = data.load_dataset()  # Load the example dataset as a dataframe
>>> est_lo=df['ci.lo']
>>> est_hi=df['ci.hi']
>>> pvalue=df['p.value']
>>> null_lo=-0.3
>>> null_hi=0.3
>>> title_lab="Leukemia Example"
>>> y_lab="Fold Change (base 10)"
>>> x_lab="Classical p-value ranking"
>>> sgpv.plot(est_lo=est_lo, est_hi=est_hi, null_lo=null_lo, null_hi=null_hi,
...            set_order=pvalue, null_pt=0, x_show=7000, outline_zone=True,
...            title_lab=title_lab, y_lab=y_lab, x_lab=x_lab )
>>> plt.yticks(ticks=np.round(np.log10(np.asarray(
...        (1/1000,1/100,1/10,1/2,1,2,10,100,1000))),2), labels=(
...                           '1/1000','1/100','1/10','1/2',1,2,10,100,1000))
>>> plt.show()
```

# Release history
* Version 1.0.3.post1: 15.07.2020:
    * Fixed a couple of formatting issues in the docstrings.
    * Cleaned the documentation of 'set_order' option of the plot function.
    * Renamed the implicit function 'power' to 'power_x' to avoid a problematic import for the risk-function. (No functional change)
* Version 1.0.3 10.07.2020:
    ## General changes 
    * Reformatted the code with autopep8 and flake8. 
    * Renamed some variables to confirm more with Python conventions. 
    * Added more descriptions based on the R-code to the documentation.  
    ## power-function
    * Fixed the display of the bonus statistic 'at 0':  Now this value is only displayed in the correct situation; the description for this value was added to the documentation.  
    ## risk-function    
    * Fixed inconsistencies/mistakes in the documentation for the risk-function.
    * Renamed the returned value of the risk-function from 'res' to 'fdcr' to reflect better the content of the variable.
    * Added a better formated output, similar to the output of the Stata version of this function.
    ## plot-function
    * Added some more input checks and added a better description of the allowed input for the option "set_order".
* Version 1.0.1 25.06.2020: Fixed incorrect imports in examples and modified code for importing the example dataset based on code found in statsmodels.datasets.utils.
* Version 1.0.0 24.06.2020: Initial release
 

# References
Blume JD, D’Agostino McGowan L, Dupont WD, Greevy RA Jr. (2018).
Second-generation p-values: Improved rigor, reproducibility, & transparency in
statistical analyses. PLoS ONE 13(3): e0188299.
[https://doi.org/10.1371/journal.pone.0188299](https://doi.org/10.1371/journal.pone.0188299)

Blume JD, Greevy RA Jr., Welty VF, Smith JR, Dupont WD (2019). An Introduction
to Second-generation p-values. The American Statistician. In press.
[https://doi.org/10.1080/00031305.2018.1537893](https://doi.org/10.1080/00031305.2018.1537893)
 

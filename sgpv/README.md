# sgpv module

This package allows to calculate Second Generation P-Values and their diagnostics in Python.
This package is a translation of the original  [sgpv R-library](https://github.com/weltybiostat/sgpv) into Python.
The same library has already been translated into [Stata](https://github.com/skbormann/stata-tools/sgpv) by me.

This module contains the following functions:
            sgpvalue   - calculate the SGPVs
            sgpower    - power functions for the SGPVs
            fdrisk    - false confirmation/discovery risks for the SGPVs
            plotsgpv   - plot the SGPVs

# Dependencies

This module depends on: 
pandas>=1.0.4
matplotlib>=3.2.1
numpy>=1.18.0
scipy>=1.3.2


These dependencies  document only under which version I tested my functions.
Older version might work as well. 


 
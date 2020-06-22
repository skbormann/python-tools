"""This module implements the calculation"""


def sgpvalue(*, null_lo, null_hi, est_lo, est_hi, inf_correction: float = 1e-5,
             warnings: bool = True):
    """
    Second-Generation p-values and delta-gaps.
    Output is still not pretty-> need to remove numpy type information
     Add user-defined expections for better error handling:
        see https://docs.python.org/3/tutorial/errors.html
    Parameters
    ----------
    * : TYPE
        DESCRIPTION.
    null_lo : array_like
        Lower bounds of the null interval(s). Values may be finite or -Inf or Inf.
        Must be of same length as null_hi.
    null_hi : array_like
        Upper bounds of the null interval(s). Values may be finite or -Inf or Inf.
        Must be of same length as null_hi.
    est_lo : array_like
        Lower bounds of interval estimates. Values may be finite or -Inf or Inf.
        Must be of same length as est_hi.
    est_hi : array_like
        Upper bounds of interval estimates. Values may be finite or -Inf or Inf.
        Must be of same length as est_lo.
    infcorrection : TYPE, optional
        A small number to denote a positive but infinitesimally small SGPV.
        Default is 1e-5. SGPVs that are infinitesimally close to 1 are assigned 1-infcorrection.
        This option can only be invoked when one of the intervals has infinite length.
    warnings : bool, optional
        Warnings toggle. Showing the warnings about potentially problematic intervals.
        Warnings are on by default.

    Returns
    -------
    pdelta : array
        second-generation p-values.
    deltagap : array
        the delta gaps, Reported as None when the corresponding
        second-generation p-value is not zero.


    Examples
    # TODO : add references to original R-code and further comments
    --------
    # Simple example for three estimated log odds ratios but the same null interval

    >>> import numpy as np
    >>> from sgpvalue import sgpvalue
    >>> lb = (np.log(1.05), np.log(1.3), np.log(0.97))
    >>> ub = (np.log(1.8), np.log(1.8), np.log(1.02))
    >>> sgpvalue(est_lo = lb, est_hi = ub, null_lo = np.log(1/1.1), null_hi = np.log(1.1))
    sgpv(pdelta=array([0.12202268, 0.        , 1.        ]), deltagap=array([None, 1.7527412602319494, None], dtype=object))
    >>> sgpvalue(est_lo = np.log(1.3), est_hi = np.inf, null_lo = np.NINF, null_hi = np.log(1.1))
    At least one interval has infinite length
    sgpv(pdelta=array([0.]), deltagap=array([0.16705408466316612], dtype=object))
    >>> sgpvalue(est_lo = np.log(1.05), est_hi = np.inf, null_lo = np.NINF, null_hi = np.log(1.1))
    At least one interval has infinite length
    sgpv(pdelta=array([0.]), deltagap=array([-0.04652001563489289], dtype=object))

    # Example t-test with simulated data
    >>> from scipy.stats import ttest_ind
    >>> from scipy.stats import norm
    >>> from scipy.stats import t
    >>> np.random.seed(1776)
    >>> x1 = norm.rvs(size=15, loc=0, scale=2)
    >>> x2 = norm.rvs(size=15, loc=3, scale=2)
    >>> se = (x1-x2).std()/np.sqrt(15)
    >>> ci1 = (x1.mean()-x2.mean()) - se*t.ppf(df=13, q=0.975)
    >>> ci2 = (x1.mean()-x2.mean()) + se*t.ppf(df=13 ,q=0.975)
    >>> sgpvalue(est_lo = ci1, est_hi = ci2, null_lo = -1, null_hi = 1)
    sgpv(pdelta=array([0.]), deltagap=array([0.41978572794859326], dtype=object))

    # Simulated two-group dichotomous data for different parameters
    >>> from scipy.stats import binom
    >>> from statsmodels.stats.proportion import proportions_ztest
    >>> np.random.seed(1492)
    >>> n = 30
    >>> p1, p2 = 0.15, 0.50
    >>> x1 = binom.rvs(1, p=p1, size=n).sum() 
    >>> x2 = binom.rvs(1, p=p2, size=n).sum()
    >>> prop1 = x1.sum()/n  # Proportion of successes
    >>> prop2 = x2.sum()/n
    >>> ci1 = (prop1 - prop2) - 1.96*np.sqrt((prop1 *(1-prop1)/n) + (prop2*(1- prop2)/n))
    >>> ci2 = (prop1 - prop2) + 1.96*np.sqrt((prop1 *(1-prop1)/n) + (prop2*(1- prop2)/n))
    >>> sgpvalue(est_lo=ci1, est_hi=ci2, null_lo=-0.2, null_hi=0.2)
    sgpv(pdelta=array([0.26748331]), deltagap=array([None], dtype=object))

    #On the log odds ratio scale
    >>> a = x1
    >>> b = x2
    >>> c = 30-x1
    >>> d = 30-x2
    >>> cior1 = np.log(a*d/(b*c)) - 1.96*np.sqrt(1/a+1/b+1/c+1/d) # Delta-method SE for log odds ratio
    >>> cior2 = np.log(a*d/(b*c)) + 1.96*np.sqrt(1/a+1/b+1/c+1/d) 
    >>> sgpvalue(est_lo=cior1, est_hi=cior2, null_lo=np.log(1/1.5), null_hi=np.log(1.5))
    sgpv(pdelta=array([0.]), deltagap=array([0.6569085789171742], dtype=object))
    """

    import numpy as np
    from termcolor import colored
    from collections import namedtuple

    # Convert inputs into np.array to emulate R behaviour
    null_lo = np.asarray(null_lo, dtype=np.float64)
    null_hi = np.asarray(null_hi, dtype=np.float64)
    est_lo = np.asarray(est_lo, dtype=np.float64)
    est_hi = np.asarray(est_hi, dtype=np.float64)

    if null_hi.size != null_lo.size:
        raise ValueError('null_lo and null_hi are of different length')

    if est_lo.size != est_hi.size:
        raise ValueError('est_lo and est_hi of different lengths')
        
    if null_lo.size != est_lo.size & null_lo.size > 1:
        raise ValueError( "Options 'null_lo' and 'null_hi' must only have one\
                 argument or exactly as many arguments as options 'est_hi' and 'est_lo'.")

    if null_lo.size == 1:
        null_lo = np.repeat(null_lo, est_lo.size)
        null_hi = np.repeat(null_hi, est_hi.size)

    # Compute Interval Lengths
    est_len = np.array(est_hi) - np.array(est_lo)
    null_len = np.array(null_hi) - np.array(null_lo)

    # Warnings -> to be added once I know how to check for these
    # -> might not be 100% correct yet
    na_any = (np.any(est_lo is None) or np.any(est_hi is None) or
              np.any(null_lo is None) or np.any(null_hi is None))

    if (na_any is True) and warnings:
        print(colored('At least one input is NA', 'red'))

    if (na_any is not True) and np.any(est_len < 0) and np.any(null_len < 0) and warnings:
        warning('At least one interval length is negative')

    if (na_any is not True) and np.any(np.isinf(abs(est_len) + abs(null_len))) and warnings:
        warning('At least one interval has infinite length')

    if (na_any is not True) and (np.any(est_len == 0) or np.any(null_len == 0)) and warnings:
        warning('At least one interval has zero length')
    # SGPV computation
    overlap = np.minimum(est_hi, null_hi) - np.maximum(est_lo, null_lo)
    overlap = np.maximum(overlap, 0)
    bottom = np.minimum(2 * null_len, est_len)
    pdelta = overlap / bottom
    # Zero-length & Infinite-length intervals
    np.where((overlap == 0), 0, pdelta)

    # Overlap finite & non-zero but bottom = Inf
    np.where(overlap != 0 & np.isfinite(overlap) & np.isinf(bottom), inf_correction, pdelta)

    # Interval estimate is a point (overlap=zero) but can be in null or equal null pt
    pdelta[(est_len == 0) & (null_len >= 0) & (est_lo >= null_lo) & (est_hi <= null_hi)] = 1

    # Null interval is a point (overlap=zero) but is in interval estimate
    pdelta[(est_len > 0) & (null_len == 0) & (est_lo <= null_lo) & (est_hi >= null_hi)] = 1/2

    # One-sided intervals with overlap; overlap == Inf & bottom==Inf
    pdelta[np.isinf(overlap) & np.isinf(bottom) & ((est_hi <= null_hi) | (est_lo >= null_lo))] = 1
    pdelta[np.isinf(overlap) & np.isinf(bottom) &
           ((est_hi > null_hi) | (est_lo < null_lo))] = 1-inf_correction

    # ## Interval estimate is entire real line and null interval is NOT entire real line
    pdelta[np.isneginf(est_lo) & np.isposinf(est_hi)] = 1/2

    # ## Null interval is entire real line
    pdelta[np.isneginf(null_lo) & np.isposinf(null_hi)] = None

    if np.any(null_lo == np.NINF) & np.any(null_hi == np.inf) and warnings:
        warning('at least one null interval is entire real line')

    # Return NA for nonsense intervals -> not working correctly yet
    pdelta[(est_lo > est_hi) | (null_lo > null_hi)] = None

    if (np.any(est_lo > est_hi) or np.any(null_lo > null_hi)) and warnings:
        warning('Some interval limits likely reversed')
    # Calculate delta gap
    deltagap = np.repeat(None, len(pdelta))
    # deltagap[! is.na(pdelta) & (pdelta == 0)] = 0
    deltagap[(pdelta is not None) & (pdelta == 0)] = 0

    gap = np.maximum(est_lo, null_lo) - np.minimum(null_hi, est_hi)

    delta = null_len / 2

    # Report unscaled delta gap if null has infinite length
    delta[null_len == np.inf] = 1

    # Report unscaled delta gap if null has length zero
    delta[null_len == 0] = 1

    dg = gap / delta

    deltagap[pdelta is not None and (pdelta == 0)] = dg[pdelta is not None and (pdelta == 0)]
    sgpv = namedtuple('sgpv', 'pdelta, deltagap')
    res = sgpv(pdelta, deltagap)
    return res

def warning(text: str):
    """
    Emulate the behavior of the R-function with the same name

    Parameters
    ----------
    text : str
        DESCRIPTION.

    Returns
    -------
    None.
    """
    from termcolor import colored
    print(colored(text, 'red'))
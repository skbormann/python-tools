def sgpvalue(*, null_lo, null_hi, est_lo, est_hi, inf_correction=1e-5, warnings=True):
    """
    sgpvalue computes Second-Generation p-values and delta-gaps
    Output is still not pretty-> need to remove numpy type information

    Parameters
    ----------
    * : TYPE
        DESCRIPTION.
    null_lo : TYPE
        DESCRIPTION.
    null_hi : TYPE
        DESCRIPTION.
    est_lo : TYPE
        DESCRIPTION.
    est_hi : TYPE
        DESCRIPTION.
    infcorrection : TYPE, optional
        DESCRIPTION. The default is 1e-5.
    warnings : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    TYPE
        DESCRIPTION.
    deltagap : TYPE
        DESCRIPTION.

    """

    #Convert inputs into np.array to emulate R behaviour
    import numpy as np
    import stop
    import math
    from termcolor import colored
    if len(null_hi) != len(null_lo):
        stop("null_lo and null_hi are of different length")
    if (len(est_lo) != len(est_hi)):
        stop('est_lo and est_hi of different lengths')

    #
    if len(null_lo) == 1:
        null_lo = np.repeat(null_lo, len(est_lo))
        null_hi = np.repeat(null_hi, len(est_hi))

    # Compute Interval Lengths
    est_len = np.array(est_hi) - np.array(est_lo)
    null_len = np.array(null_hi) - np.array(null_lo)

    # Warnings -> to be added once I know how to check for these
    na_any = (np.any(est_lo is None) or np.any(est_hi is None) or 
              np.any(null_lo is None) or np.any(null_hi is None))

    if ((na_any==True) and warnings):
        print(colored('At least one input is NA', 'red'))

    if (na_any==False) and np.any(est_len<0) and np.any(null_len<0) and warnings :
        warning('At least one interval length is negative')

    if ((na_any==False) and np.any((abs(est_len)+abs(null_len))is np.inf) and warnings) :
        warning('At least one interval has infinite length') 

    if ((na_any==False) and np.any(est_len==0,null_len==0) and warnings) :
        warning('At least one interval has zero length') 
    # SGPV computation
    overlap = np.minimum(est_hi, null_hi) - np.maximum(est_lo, null_lo)
    overlap = np.maximum(overlap, 0)

    bottom = np.minimum(2 * null_len, est_len)

    pdelta = overlap / bottom

    # Zero-length & Infinite-length intervals -> to be added once I know how to check for these
    ## Overwrite NA and NaN due to bottom = Inf
    pdelta[overlap==0] = 0

    ## Overlap finite & non-zero but bottom = Inf
    pdelta[overlap!=0 and math.isfinite(overlap) and np.isinf(bottom)] = inf_correction

    ## Interval estimate is a point (overlap=zero) but can be in null or equal null pt
    pdelta[est_len==0  and null.len>=0  and est_lo>=null_lo and est_hi<=null_hi] = 1

    ## Null interval is a point (overlap=zero) but is in interval estimate
    pdelta[est_len>0  and null.len==0  and est_lo<=null_lo  and est_hi>=null_hi] = 1/2

    ## One-sided intervals with overlap; overlap == Inf & bottom==Inf
    pdelta[np.isinf(overlap) and np.isinf(bottom) and ((est_hi<=null_hi) or (est_lo>=null_lo))] = 1
    pdelta[np.isinf(overlap) and np.isinf(bottom) and ((est_hi>null_hi) or (est_lo<null_lo))] = 1-inf_correction

    ## Interval estimate is entire real line and null interval is NOT entire real line
    pdelta[est_lo==np.NINF  and est_hi==np.Inf] = 1/2

    ## Null interval is entire real line
    pdelta[null_lo==np.NINF and null_hi==np.inf] = None

    if (np.any(null_lo==np.NINF) and np.any(null_hi==np.inf) and warnings): 
        warning('at least one null interval is entire real line') 

    ## Return NA for nonsense intervals
    pdelta[(est_lo>est_hi) or (null_lo>null_hi)] = None

    if (np.any(est_lo>est_hi) or np.any(null_lo>null_hi)) and warnings:
        warning('Some interval limits likely reversed')
    # Calculate delta gap
    deltagap = np.repeat(None, len(pdelta))
    #deltagap[! is.na(pdelta) & (pdelta == 0)] = 0
    deltagap[(pdelta is not None) & (pdelta == 0)] = 0

    gap = np.maximum(est_lo, null_lo) - np.minimum(null_hi, est_hi)

    delta = null_len / 2

    # Report unscaled delta gap if null has infinite length
    delta[null_len == np.inf] = 1

    # Report unscaled delta gap if null has length zero
    delta[null_len == 0] = 1

    dg = gap / delta

    #deltagap[! is.None(pdelta) & (pdelta == 0)] = dg[! is.na(pdelta) & (pdelta == 0)]
    deltagap[pdelta is not None and (pdelta == 0)] = dg[pdelta is not None and (pdelta == 0)]

    return {'pdelta':pdelta, 'deltagap':deltagap}

def warning(text):
    """
    Emulate the behavior of the R-function with the same name

    Parameters
    ----------
    text : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    from termcolor import colored
    print(colored(text,'red'))
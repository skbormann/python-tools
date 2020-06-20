# Add type hints for parameters
def sgpower(*, true: float, null_lo, null_hi, std_err: float = 1,
            interval_type: str, interval_level):
    """
        Compute power/type I error for Second-Generation p-values approach.
    For now only single input values are fully supported
    vector inputs are still error-prone
    No error checks
    Add user-defined expections for better error handling:
        see https://docs.python.org/3/tutorial/errors.html

    Parameters
    ----------
    * : TYPE
        DESCRIPTION.
    true : float
        DESCRIPTION.
    null_lo : TYPE
        DESCRIPTION.
    null_hi : TYPE
        DESCRIPTION.
    std_err : float, optional
        DESCRIPTION. The default is 1.
    interval_type : str
        DESCRIPTION.
    interval_level : TYPE
        DESCRIPTION.

    Raises
    ------
    to
        DESCRIPTION.

    Returns
    -------
    res : TYPE
        DESCRIPTION.

    Examples
    --------
    sgpower(true=2, null_lo-1, null_hi=1. interval_type='confidence', interval_level=0.05)

    """
    import numpy as np
    from scipy import integrate
    from scipy.stats import norm
    from collections import namedtuple # for better naming of the output
    from stop import stop
    #Need to learn error handling
    if interval_type not in ['confidence', 'likelihood']:
        stop("Parameter `intervaltype` must be one of the following: \
             \n   'confidence' \n  'likelihood' \n  \
             \ (credible not currently supported for power calculations)")

    if interval_type == 'confidence':
        Z = norm.ppf(1-interval_level/2)
    if interval_type == 'likelihood':
        Z = norm.ppf(1-2*norm.cdf(-np.sqrt(2*np.log(1/interval_level)))/2)

    ## P(SGPV=0 | true )    (see Blume et al. (2018) eq.(S4) for CI/LSI)
    power0 = norm.cdf(null_lo/std_err - true/std_err - Z)
    + norm.cdf(-null_hi/std_err + true/std_err - Z)

    ## P(SGPV=1 | true )    (see Blume et al. (2018) eq.(S7) for CI/LSI)
    # -> only for symmetric null hypothesis
    if (null_hi-null_lo) >= 2*Z*std_err:
        power1 = norm.cdf(null_hi/std_err - true/std_err - Z)
        - norm.cdf(null_lo/std_err - true/std_err + Z)

    if (null_hi-null_lo) < 2*Z*std_err:
        power1 = 0

    ## P(0<SGPV<1 | true)   (see Blume et al. (2018) eq.(S8, S9) for CI/LSI)
    # -> only for symmetric null hypothesis
    if (null_hi-null_lo) <= 2*Z*std_err:
        powerinc = 1 - norm.cdf(null_lo/std_err - true/std_err - Z)
        - norm.cdf(-null_hi/std_err + true/std_err - Z)

    if (null_hi-null_lo) > 2*Z*std_err:
        powerinc = 1 - (norm.cdf(null_lo/std_err - true/std_err - Z)
                        + norm.cdf(-null_hi/std_err + true/std_err - Z))
        - (norm.cdf(null_hi/std_err - true/std_err - Z)
           - norm.cdf(null_lo/std_err - true/std_err + Z))

    ## check -> works only for scalar input
    #Need to find out which exception to raise to emulate behaviour of R's warning()-function
    #Need to find out how the index
    if round(power0 + powerinc + power1, 7) != 1:
        print('error: power0+powerinc+power1 != 1 for indices ',
              round(power0 + powerinc + power1, 7) != 1)
    ## bonus: type I error summaries
    pow0 = lambda x: norm.cdf(null_lo/std_err - x/std_err - Z) + norm.cdf(
                            -null_hi/std_err + x/std_err - Z)
    minI = pow0((null_lo + null_hi)/2)
    maxI = pow0(null_lo)
    avgI = 1/(null_hi - null_lo)*integrate.quad(pow0, null_lo, null_hi)[0]
    typeI = (minI, maxI, avgI)
    if null_lo <= 0 & 0 <= null_hi:
        #Not quite the output I want
        #typeI = ('at 0 =', round(pow0(0),7), 'min =', round(minI,7), 'max =',round(maxI,7), 'mean =', round(avgI,7))
        TypeI = namedtuple('typeI','at0, min, max, mean')
        typeI = TypeI(round(pow0(0), 7), round(minI, 7), round(maxI, 7), round(avgI, 7))
    ## more bonus: P(inc | null) (basically analogous to type II error but for confirmation)
    ## (TBD)
    #Displaying of bonus statistics not correct yet
   # print('poweralt' = power0, 'powerinc' = powerinc, 'powernull' = power1, 'type I error summaries' = typeI)
    print('poweralt =', round(power0, 7), 'powerinc =', round(powerinc, 7), \
          'powernull = ', round(power1, 7), '\n type I error summaries: \n', \
              'at 0 =', round(pow0(0), 7), 'min =', round(minI,7),  \
                  'max =', round(maxI,7), 'mean =', round(avgI,7))
    sgpow = namedtuple('sgpower', ['poweralt', 'powerinc','powernull','typeI'])
    res = sgpow(round(power0, 7), round(powerinc, 7), round(power1, 7), typeI)
    return res

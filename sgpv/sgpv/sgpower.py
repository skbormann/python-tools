"""Compute power/type I error for Second-Generation p-values approach."""


def sgpower(*, true: float, null_lo, null_hi, std_err: float = 1,
            interval_type: str, interval_level: float):
    """
        Compute power/type I error for Second-Generation p-values approach.
    For now only single input values are fully supported
    vector inputs are still error-prone
    No error checks
    Add user-defined expections for better error handling:
        see https://docs.python.org/3/tutorial/errors.html

    Parameters
    ----------
    true : float
        The true value for the parameter of interest at which to calculate power. 
        This is on the absolute scale of the parameter, and not the standard
        deviation or standard error scale.
    null_lo : TYPE
        The lower bound of the indifference zone (null interval) upon which the 
        second-generation p-value is based.
    null_hi : TYPE
        The upper bound of the indifference zone (null interval) upon which the
        second-generation p-value is based.
    std_err : float, optional
        Standard error for the distribution of the estimator for the parameter of interest. 
        This is the standard deviation for the estimator, not the standard 
        deviation parameter for the data itself. This will be a function of the 
        sample size(s). The default is 1.
    interval_type : str
        Class of interval estimate used. This determines the functional form of
        the power function. Options are 'confidence' for a (1-$\alpha$)100\%
        confidence interval and 'likelihood' for a 1/k likelihood support interval.
    interval_level : float
        Level of interval estimate. If inttype is 'confidence', the level is $\alpha$. 
        If interval_type is 'likelihood', the level is 1/k (not k).

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
    >>> from sgpower import sgpower
    >>> sgpower(true=2, null_lo=-1, null_hi=1, std_err = 1,
    ...        interval_type='confidence', interval_level=0.05)
    poweralt = 0.168537 powerinc = 0.831463 powernull =  0
    type I error summaries:
    at 0 = 0.0030768 min = 0.0030768 max = 0.0250375 mean = 0.0094374
    >>> sgpower(true=0, null_lo=-1, null_hi=1, std_err = 1,
    ...         interval_type='confidence', interval_level=0.05)
    poweralt = 0.0030768 powerinc = 0.9969232 powernull =  0
    type I error summaries:
    at 0 = 0.0030768 min = 0.0030768 max = 0.0250375 mean = 0.0094374
    
    #Plot the power curve
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> sigma, n = 5, 20
    >>> se = sigma/np.sqrt(n)
    >>> theta = np.arange(-10,10,0.1)
    >>> fig, ax = plt.subplots()
    >>> power = []
    >>> for i in range(len(theta)):
    ...     power.append( sgpower(true=theta[i], null_lo=-1, null_hi=1, std_err=se,
    ...                interval_type='confidence', interval_level=0.05).poweralt)
    >>> ax.plot(theta, power)
    >>> ax.set_xlabel('theta')
    >>> ax.set_ylabel('power')
    """
    import numpy as np
    from scipy import integrate
    from scipy.stats import norm
    from collections import namedtuple  # for better naming of the output
    # Need to learn error handling
    if interval_type not in ['confidence', 'likelihood']:
        raise ValueError("Parameter 'interval_type' must be one of the following:\
             \n   'confidence'  or  'likelihood'")

    if interval_type == 'confidence':
        Z = norm.ppf(1-interval_level/2)
    if interval_type == 'likelihood':
        Z = norm.ppf(1-2*norm.cdf(-np.sqrt(2*np.log(1/interval_level)))/2)

    # P(SGPV=0 | true )    (see Blume et al. (2018) eq.(S4) for CI/LSI)
    power0 = norm.cdf(null_lo/std_err - true/std_err - Z) + norm.cdf(
        -null_hi/std_err + true/std_err - Z)

    # P(SGPV=1 | true )    (see Blume et al. (2018) eq.(S7) for CI/LSI)
    # -> only for symmetric null hypothesis
    if (null_hi-null_lo) >= 2*Z*std_err:
        power1 = norm.cdf(null_hi/std_err - true/std_err - Z) - norm.cdf(
            null_lo/std_err - true/std_err + Z)

    if (null_hi-null_lo) < 2*Z*std_err:
        power1 = 0

    # P(0<SGPV<1 | true)   (see Blume et al. (2018) eq.(S8, S9) for CI/LSI)
    # -> only for symmetric null hypothesis
    if (null_hi-null_lo) <= 2*Z*std_err:
        powerinc = 1 - norm.cdf(null_lo/std_err - true/std_err - Z) - norm.cdf(
            -null_hi/std_err + true/std_err - Z)

    if (null_hi-null_lo) > 2*Z*std_err:
        powerinc = 1 - (norm.cdf(null_lo/std_err - true/std_err - Z)
                        + norm.cdf(-null_hi/std_err + true/std_err - Z))-(
                        norm.cdf(null_hi/std_err - true/std_err - Z) - norm.cdf(
                            null_lo/std_err - true/std_err + Z))

    # check -> works only for scalar input
    # Need to find out which exception to raise to emulate behaviour of R's warning()-function
    # Need to find out how to the index
    if round(power0 + powerinc + power1, 2) != 1:
        print('error: power0+powerinc+power1 != 1 for indices ',
              round(power0 + powerinc + power1, 2) != 1)
    # bonus: type I error summaries
    pow0 = lambda x: norm.cdf(null_lo/std_err - x/std_err - Z) + norm.cdf(
        -null_hi/std_err + x/std_err - Z)
    minI = pow0((null_lo + null_hi)/2)
    maxI = pow0(null_lo)
    avgI = 1/(null_hi - null_lo)*integrate.quad(pow0, null_lo, null_hi)[0]
    typeI = (minI, maxI, avgI)
    if null_lo <= 0 & 0 <= null_hi:
        # Not quite the output I want
        # typeI = ('at 0 =', round(pow0(0),7), 'min =', round(minI,7), 'max =',round(maxI,7), 'mean =', round(avgI,7))
        TypeI = namedtuple('typeI', 'at0, min, max, mean')
        typeI = TypeI(round(pow0(0), 7), round(minI, 7), round(maxI, 7), round(avgI, 7))
    # more bonus: P(inc | null) (basically analogous to type II error but for confirmation)
    # (TBD)
    # Displaying of bonus statistics not correct yet
    print('poweralt =', round(power0, 7), 'powerinc =', round(powerinc, 7),
          'powernull = ', round(power1, 7), '\n type I error summaries: \n',
          'at 0 =', round(pow0(0), 7), 'min =', round(minI, 7),
          'max =', round(maxI, 7), 'mean =', round(avgI, 7))
    sgpow = namedtuple('sgpower', ['poweralt', 'powerinc', 'powernull', 'typeI'])
    res = sgpow(round(power0, 7), round(powerinc, 7), round(power1, 7), typeI)
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
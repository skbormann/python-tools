"""Implements the calculation of Second Generation P-Values and their associated diagnostics.

Author: Sven-Kristjan Bormann
This module implements the calculation of Second Generation P-Values and
their associated diagnostics.

This module contains the following functions:

            value    - calculate the SGPVs.

            power    - power functions for the SGPVs.

            risk     - false confirmation/discovery risks for the SGPVs.

            plot     - Plot interval estimates according to SGPV rankings.

The sgpv-module is a translation of the same named R-library of the same name.
The original R-library can be found here https://github.com/weltybiostat/sgpv

References:
Blume JD, D’Agostino McGowan L, Dupont WD, Greevy RA Jr. (2018).
Second-generation p-values: Improved rigor, reproducibility, & transparency in
statistical analyses. PLoS ONE 13(3): e0188299.
https://doi.org/10.1371/journal.pone.0188299

Blume JD, Greevy RA Jr., Welty VF, Smith JR, Dupont WD (2019). An Introduction
to Second-generation p-values. The American Statistician. In press.
https://doi.org/10.1080/00031305.2018.1537893
"""


def value(*, null_lo, null_hi, est_lo, est_hi, inf_correction: float = 1e-5,
          warnings: bool = True):
    """Second-Generation p-values and delta-gaps.

    This function computes the second-generation p-value (SGPV) and
    its associated delta gaps, as introduced in Blume et al. (2018).

    Parameters
    ----------
    null_lo : array_like
        Lower bounds of the null interval(s). Values may be finite or -Inf or Inf.
        Must have the same number of elements as null_hi.
    null_hi : array_like
        Upper bounds of the null interval(s). Values may be finite or -Inf or Inf.
        Must have the same number of elements as null_hi.
    est_lo : array_like
        Lower bounds of interval estimates. Values may be finite or -Inf or Inf.
        Must have the same number of elements as est_hi.
    est_hi : array_like
        Upper bounds of interval estimates. Values may be finite or -Inf or Inf.
        Must have the same number of elements as est_lo.
    infcorrection : float, optional
        A small number to denote a positive but infinitesimally small SGPV.
        Default is 1e-5. SGPVs that are infinitesimally close to 1
        are assigned 1-infcorrection.
        This option can only be invoked when one of the intervals
        has infinite length.
    warnings : bool, optional
        Warnings toggle. Showing the warnings about potentially
        problematic intervals.
        Warnings are on by default.

    Raises
    ------
    ValueError
        Indicates that some value was outside of the expected range or
        outside of the accepted options.

    Returns
    -------
    pdelta : numpy_array
        Second-generation p-values.
    deltagap : numpy_array
        The delta gaps, Reported as None when the corresponding
        second-generation p-value is not zero.

    Examples
    --------
    Simple example for three estimated log odds ratios but the same null interval.

    >>> import numpy as np
    >>> from sgpv import sgpv
    >>> lb = (np.log(1.05), np.log(1.3), np.log(0.97))
    >>> ub = (np.log(1.8), np.log(1.8), np.log(1.02))
    >>> sgpv.value(est_lo = lb, est_hi = ub,
    ...             null_lo = np.log(1/1.1), null_hi = np.log(1.1))
    sgpv(pdelta=array([0.1220227, 0.        , 1.        ]),
    ...     deltagap=array([None, 1.7527413, None], dtype=object))
    >>> sgpv.value(est_lo = np.log(1.3), est_hi = np.inf,
    ...             null_lo = np.NINF, null_hi = np.log(1.1))
    At least one interval has infinite length
    sgpv(pdelta=array([0.]), deltagap=array([0.1670541], dtype=object))
    >>> sgpv.value(est_lo = np.log(1.05), est_hi = np.inf,
    ...             null_lo = np.NINF, null_hi = np.log(1.1))
    At least one interval has infinite length
    sgpv(pdelta=array([0.]), deltagap=array([-0.04652], dtype=object))

    Example t-test with simulated data

    >>> from scipy.stats import ttest_ind
    >>> from scipy.stats import norm
    >>> from scipy.stats import t
    >>> np.random.seed(1776)
    >>> x1 = norm.rvs(size=15, loc=0, scale=2)
    >>> x2 = norm.rvs(size=15, loc=3, scale=2)
    >>> se = (x1-x2).std()/np.sqrt(15)
    >>> ci1 = (x1.mean()-x2.mean()) - se*t.ppf(df=13, q=0.975)
    >>> ci2 = (x1.mean()-x2.mean()) + se*t.ppf(df=13 ,q=0.975)
    >>> sgpv.value(est_lo = ci1, est_hi = ci2,
    ...               null_lo = -1, null_hi = 1)
    sgpv(pdelta=array([0.]), deltagap=array([0.3000322], dtype=object))

    Simulated two-group dichotomous data for different parameters

    >>> from scipy.stats import binom
    >>> np.random.seed(1492)
    >>> n = 30
    >>> p1, p2 = 0.15, 0.50
    >>> x1 = binom.rvs(1, p=p1, size=n).sum()
    >>> x2 = binom.rvs(1, p=p2, size=n).sum()
    >>> prop1 = x1.sum()/n  # Proportion of successes
    >>> prop2 = x2.sum()/n
    >>> ci1 = (prop1 - prop2) - 1.96*np.sqrt((prop1 *(1-prop1)/n)
    ...                                    + (prop2*(1- prop2)/n))
    >>> ci2 = (prop1 - prop2) + 1.96*np.sqrt((prop1 *(1-prop1)/n)
    ...                                    + (prop2*(1- prop2)/n))
    >>> sgpv.value(est_lo=ci1, est_hi=ci2, null_lo=-0.2, null_hi=0.2)
    sgpv(pdelta=array([0.2756205]),
         deltagap=array([None], dtype=object))

    On the log odds ratio scale

    >>> a = x1
    >>> b = x2
    >>> c = 30-x1
    >>> d = 30-x2
    >>> cior1 = np.log(a*d/(b*c)) - 1.96*np.sqrt(1/a+1/b+1/c+1/d)
    >>> cior2 = np.log(a*d/(b*c)) + 1.96*np.sqrt(1/a+1/b+1/c+1/d)
    >>> sgpv.value(est_lo=cior1, est_hi=cior2,
    ...          null_lo=np.log(1/1.5), null_hi=np.log(1.5))
    sgpv(pdelta=array([0.]), deltagap=array([0.65691], dtype=object))
    """
    from collections import namedtuple

    import numpy as np

    # Convert inputs into np.array to emulate R behaviour.
    null_lo = np.asarray(null_lo, dtype=np.float64)
    null_hi = np.asarray(null_hi, dtype=np.float64)
    est_lo = np.asarray(est_lo, dtype=np.float64)
    est_hi = np.asarray(est_hi, dtype=np.float64)

    if null_hi.size != null_lo.size:
        raise ValueError('null_lo and null_hi are of different lengths.')

    if est_lo.size != est_hi.size:
        raise ValueError('est_lo and est_hi are of different lengths.')

    if null_lo.size != est_lo.size & null_lo.size > 1:
        raise ValueError("'null_lo' and 'null_hi' must only have one argument\
                     or exactly as many arguments as 'est_hi' and 'est_lo'.")

    if null_lo.size == 1:
        null_lo = np.repeat(null_lo, est_lo.size)
        null_hi = np.repeat(null_hi, est_hi.size)

    # Compute interval lengths.
    est_len = np.array(est_hi) - np.array(est_lo)
    null_len = np.array(null_hi) - np.array(null_lo)

    # Warnings -> might not be 100% correct yet
    na_any = (np.any(est_lo is None) or np.any(est_hi is None) or
              np.any(null_lo is None) or np.any(null_hi is None))

    if (na_any is True) and warnings:
        print('At least one input is None.')

    if (na_any is not True) and np.any(
            est_len < 0) and np.any(null_len < 0) and warnings:
        print('At least one interval length is negative.')

    if (na_any is not True) and np.any(
            np.isinf(abs(est_len) + abs(null_len))) and warnings:
        print('At least one interval has infinite length.')

    if (na_any is not True) and (np.any(est_len == 0)
                                 or np.any(null_len == 0)) and warnings:
        print('At least one interval has zero length.')

    # SGPV computation
    overlap = np.minimum(est_hi, null_hi) - np.maximum(est_lo, null_lo)
    overlap = np.maximum(overlap, 0)
    bottom = np.minimum(2 * null_len, est_len)
    pdelta = np.round(overlap / bottom, 7)

    # Zero-length & Infinite-length intervals
    np.where((overlap == 0), 0, pdelta)

    # Overlap finite & non-zero but bottom = Inf
    np.where(overlap != 0 & np.isfinite(overlap) &
             np.isinf(bottom), inf_correction, pdelta)

    # Interval estimate is a point (overlap=zero) but can be in null or equal
    # null pt
    pdelta[(est_len == 0) & (null_len >= 0) & (
        est_lo >= null_lo) & (est_hi <= null_hi)] = 1

    # Null interval is a point (overlap=zero) but is in interval estimate
    pdelta[(est_len > 0) & (null_len == 0) & (
        est_lo <= null_lo) & (est_hi >= null_hi)] = 1 / 2

    # One-sided intervals with overlap; overlap == Inf & bottom==Inf
    pdelta[np.isinf(overlap) & np.isinf(bottom) & (
        (est_hi <= null_hi) | (est_lo >= null_lo))] = 1
    pdelta[np.isinf(overlap) & np.isinf(bottom) &
           ((est_hi > null_hi) | (est_lo < null_lo))] = 1 - inf_correction

    # Interval estimate is entire real line and null interval is NOT entire
    # real line
    pdelta[np.isneginf(est_lo) & np.isposinf(est_hi)] = 1 / 2

    # Null interval is entire real line
    pdelta[np.isneginf(null_lo) & np.isposinf(null_hi)] = None

    if np.any(null_lo == np.NINF) & np.any(null_hi == np.inf) and warnings:
        print('At least one null interval is entire real line.')

    # Return None for nonsense intervals -> not working correctly yet
    pdelta[(est_lo > est_hi) | (null_lo > null_hi)] = None

    if (np.any(est_lo > est_hi) or np.any(null_lo > null_hi)) and warnings:
        print('Some interval limits likely reversed.')

    # Calculate delta gap
    deltagap = np.repeat(None, len(pdelta))
    deltagap[(pdelta is not None) & (pdelta == 0)] = 0

    gap = np.maximum(est_lo, null_lo) - np.minimum(null_hi, est_hi)

    delta = null_len / 2

    # Report unscaled delta gap if null has infinite length
    delta[null_len == np.inf] = 1

    # Report unscaled delta gap if null has length zero
    delta[null_len == 0] = 1

    dg = np.round(gap / delta, 7)

    deltagap[pdelta is not None and (pdelta == 0)] = dg[
        pdelta is not None and (pdelta == 0)]
    sgpv = namedtuple('sgpv', 'pdelta, deltagap')
    return sgpv(pdelta, deltagap)


def power(*, true: float, null_lo, null_hi, std_err: float = 1,
          interval_type: str, interval_level: float, no_print: bool = False):
    """Compute power/type I error for Second-Generation p-values approach.

    Calculate power and type I error values from significance testing based
    on second-generation p-values as the inferential metric.

    Note
    ----
    For now only single input values are fully supported.

    Parameters
    ----------
    true : float
        The true value for the parameter of interest at which to calculate
        the power function.
        This is on the absolute scale of the parameter, and not the standard
        deviation or standard error scale.
    null_lo : float
        The lower bound of the indifference zone (null interval) upon which the
        second-generation p-value is based.
    null_hi : float
        The upper bound of the indifference zone (null interval) upon which the
        second-generation p-value is based.
    std_err : float, optional
        Standard error for the distribution of the estimator for the parameter
        of interest. This is the standard deviation for the estimator,
        not the standard deviation parameter for the data itself. This will be
        a function of the sample size(s).
        The default is 1.
    interval_type : str
        Type of interval estimate used. This determines the functional form of
        the power function. Options are 'confidence' for a (1-α)100% confidence
        interval and 'likelihood' for a 1/k likelihood support interval.
    interval_level : float
        Level of interval estimate. If inttype is 'confidence', the level is α.
        If interval_type is 'likelihood', the level is 1/k (not k).
    no_print : bool, optional
        Disable printing of a better formatted output. Default is False,
        so that many lines maybe printed if the power-function is used repeatedly.

    Raises
    ------
    ValueError
        Indicates that some value was outside of the expected range or
        outside of the accepted options.

    Returns
    -------
    sgpow : tuple
        A list containing the following components:
        poweralt:
            Probability of SGPV = 0 calculated assuming the parameter
                  is equal to true. That is, poweralt = P(SGPV = 0 | θ = true).
        powerinc:
            Probability of 0 < SGPV < 1 calculated assuming the parameter
                  is equal to true. That is, poweralt = P(0 < SGPV < 1 | θ = true).
        powernull:
            Probability of SGPV = 1 calculated assuming the parameter
                  is equal to true. That is, poweralt = P(SGPV = 1 | θ = true).

        type I error summaries:
                                Named vector that includes different ways
                                the type I error may be summarized for an
                                interval null hypothesis.
                    min:
                        min is the minimum type I error over the range
                        (null_lo, null_hi), which occurs at the midpoint of
                        (null_lo, null_hi).
                    max:
                        is the maximum type I error over the range
                        (null_lo, null_hi), which occurs at the boundaries of
                        the null hypothesis, null_lo and null_hi.
                    mean:
                        is the average type I error (unweighted) over
                        the range (null_lo, null_hi).

                        If 0 is included in the null hypothesis region,
                        then `type I error summaries` also contains 'at 0'
                    pow0:
                        is the type I error calculated assuming the true parameter
                        value θ is equal to 0.

    Examples
    --------
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

    Plot the power curve

    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> sigma, n = 5, 20
    >>> se = sigma/np.sqrt(n)
    >>> theta = np.arange(-10,10,0.1)
    >>> fig, ax = plt.subplots()
    >>> power = []
    >>> for i in range(len(theta)):
    ...     power.append(
                sgpv.power(true=theta[i], null_lo=-1, null_hi=1,
    ...                 std_err=se, interval_type='confidence',
    ...                 interval_level=0.05, no_print=True ).poweralt)
    >>> ax.plot(theta, power)
    >>> ax.set_xlabel('theta')
    >>> ax.set_ylabel('power')
    >>> plt.show()
    """
    from collections import namedtuple

    import numpy as np
    from scipy import integrate
    from scipy.stats import norm

    # Input checks
    if interval_type not in ['confidence', 'likelihood']:
        raise ValueError("'interval_type' must be one of the following:\
             \n   'confidence'  or  'likelihood'")

    inputs = ['null_lo', 'null_hi', 'true', 'std_err', 'interval_level']
    for index, val in enumerate(
            [null_lo, null_hi, true, std_err, interval_level]):
        if np.asarray(val).size != 1:
            raise ValueError(f"The argument {inputs[index]} must be a scalar.")

    if interval_type == 'confidence':
        Z = norm.ppf(1 - interval_level / 2)
    if interval_type == 'likelihood':
        Z = norm.ppf(
            1 - 2 * norm.cdf(-np.sqrt(2 * np.log(1 / interval_level))) / 2)

    # P(SGPV=0 | true )    (see Blume et al. (2018) eq.(S4) for CI/LSI)
    power0 = norm.cdf(null_lo / std_err - true / std_err - Z) + norm.cdf(
        -null_hi / std_err + true / std_err - Z)

    # P(SGPV=1 | true )    (see Blume et al. (2018) eq.(S7) for CI/LSI)
    # -> only for symmetric null hypothesis
    if (null_hi - null_lo) >= 2 * Z * std_err:
        power1 = norm.cdf(null_hi / std_err - true / std_err - Z) - norm.cdf(
            null_lo / std_err - true / std_err + Z)

    if (null_hi - null_lo) < 2 * Z * std_err:
        power1 = 0

    # P(0<SGPV<1 | true)   (see Blume et al. (2018) eq.(S8, S9) for CI/LSI)
    # -> only for symmetric null hypothesis
    if (null_hi - null_lo) <= 2 * Z * std_err:
        powerinc = 1 - norm.cdf(null_lo / std_err - true / std_err - Z) - \
            norm.cdf(-null_hi / std_err + true / std_err - Z)

    if (null_hi - null_lo) > 2 * Z * std_err:
        powerinc = 1 - (norm.cdf(null_lo / std_err - true / std_err - Z)
                        + norm.cdf(-null_hi / std_err + true / std_err - Z)) - (
            norm.cdf(null_hi / std_err - true / std_err - Z) - norm.cdf(
                null_lo / std_err - true / std_err + Z))

    # check -> works only for scalar input
    # Need to find out how to get the index
    if round(power0 + powerinc + power1, 2) != 1:
        print('error: power0+powerinc+power1 != 1',
              round(power0 + powerinc + power1, 2) != 1)

    # bonus: type I error summaries
    def pow0(x): return norm.cdf(null_lo / std_err - x / std_err -\
             Z) + norm.cdf(-null_hi / std_err + x / std_err - Z)
    minI = pow0((null_lo + null_hi) / 2)
    maxI = pow0(null_lo)
    avgI = 1 / (null_hi - null_lo) * integrate.quad(pow0, null_lo, null_hi)[0]

    if not null_lo <= 0 & 0 <= null_hi:
        TypeI = namedtuple('typeI', 'min, max, mean')
        typeI = TypeI(round(minI, 7), round(maxI, 7), round(avgI, 7))

    if null_lo <= 0 & 0 <= null_hi:
        TypeI = namedtuple('typeI', 'at0, min, max, mean')
        typeI = TypeI(round(pow0(0), 7), round(minI, 7),
                      round(maxI, 7), round(avgI, 7))

    # Print bonus statistics only if requested
    if no_print is False:
        if null_lo <= 0 & 0 <= null_hi:
            print('poweralt =', round(power0, 7),
                  'powerinc =', round(powerinc, 7),
                  'powernull = ', round(power1, 7),
                  '\n type I error summaries: \n', 'at 0 =', round(pow0(0), 7),
                  'min =', round(minI, 7),
                  'max =', round(maxI, 7),
                  'mean =', round(avgI, 7))
        if not null_lo <= 0 & 0 <= null_hi:
            print('poweralt =', round(power0, 7),
                  'powerinc =', round(powerinc, 7),
                  'powernull = ', round(power1, 7),
                  '\n type I error summaries: \n',
                  'min =', round(minI, 7),
                  'max =', round(maxI, 7),
                  'mean =', round(avgI, 7))
    sgpow = namedtuple(
        'sgpower', [
            'poweralt', 'powerinc', 'powernull', 'typeI'])
    return sgpow(round(power0, 7), round(powerinc, 7), round(power1, 7), typeI)


def risk(
        *,
        sgpval: int = 0,
        null_lo: float,
        null_hi: float,
        std_err: float,
        interval_type: str,
        interval_level: float,
        null_weights: str,
        null_space,
        alt_weights: str,
        alt_space,
        pi0: float = 0.5):
    """Compute False Discovery or Confirmation risk for the Second-Generation p-values approach.

    This function computes the false discovery risk (sometimes called the
    "empirical bayes FDR") for a second-generation p-value of 0, or
    the false confirmation risk for a second-generation p-value of 1.

    Parameters
    ----------
    sgpval : int, optional
        The observed second-generation p-value.
        Default is 0, which gives the false discovery risk.
        Setting it to 1 gives the false confirmation risk.
    null_lo : float
        The lower bound of the indifference zone (null interval) upon which the
        second-generation p-value was based.
    null_hi : float
        The upper bound of the indifference zone (null interval) upon which the
        second-generation p-value was based.
    std_err : float
        Standard error of the point estimate.
    interval_type : str
        Type of interval estimate used.
        This determines the functional form of the power function.
        Options are 'confidence' for a (1-α)100 % confidence interval and
        'likelihood' for a 1/k likelihood support interval.
    interval_level : float
        Level of interval estimate. If inttype is 'confidence', the level is α.
        If interval_type is 'likelihood', the level is 1/k (not k).
    null_weights : str
        Probability distribution for the null parameter space.
        Options are 'Point', 'Uniform', and 'TruncNormal'.
    null_space : array_like
        Support of the null probability distribution.
        If null_weights is 'Point', then null_space is a scalar.
        If null_weights is 'Uniform' or 'TruncNormal', then
        null_space is an array or list of length two.
    alt_weights : str
        Probability distribution for the alternative parameter space.
        Options are 'Point', 'Uniform', and 'TruncNormal'.
    alt_space : array_like
        Support of the alternative probability distribution.
        If alt_weights is 'Point', then alt_space is a scalar.
        If alt_weights is 'Uniform' or 'TruncNormal',
        then alt_space is an array or list of length two.
    pi0 : float, optional
        Prior probability of the null hypothesis. The default is 0.5.
        This value can be only between 0 and 1 (exclusive).
        A prior probability outside of this interval is not sensible.
        The default value assumes that both hypotheses are equally likely.

    Details
    -------
    When possible, one should compute the second-generation p-value and
    FDR/FCR on a scale that is symmetric about the null hypothesis.
    For example, if the parameter of interest is an odds ratio,
    inputs std_err, null_lo, null_hi, null_space, and alt_space are
    typically on the log scale.

    If TruncNormal is used for null_weights,
    then the distribution used is a truncated Normal distribution with
    mean equal to the midpoint of null_space, and standard deviation equal to
    std_err, truncated to the support of null_space.
    If TruncNormal is used for alt_weights,
    then the distribution used is a truncated Normal distribution with
    mean equal to the midpoint of alt_space, and standard deviation equal to
    std_err, truncated to the support of alt.space.
    Further customization of these parameters for the truncated Normal are
    not possible.

    Raises
    ------
    ValueError
        Indicates that some value was outside of the expected range or
        outside of the accepted options.


    Returns
    -------
    fdcr: float
        Numeric scalar representing the False discovery risk (FDR) or
        false confirmation risk (FCR) for the observed second-generation p-value.

        If sgpval = 0, the function returns false discovery risk (FDR).
        If sgpval = 1, the function returns false confirmation risk (FCR).


    Examples
    --------
    >>> from sgpv import sgpv
    >>> import numpy as np
    >>> from scipy.stats import norm

    False discovery risk with 95% confidence level

    >>> sgpv.risk(sgpval = 0, null_lo = np.log(1/1.1),
                  null_hi = np.log(1.1), std_err = 0.8,
                  null_weights = 'Uniform',
               null_space = (np.log(1/1.1), np.log(1.1)),
               alt_weights = 'Uniform',
               alt_space = (2 + 1*norm.ppf(1-0.05/2)*0.8,
                            2 - 1*norm.ppf(1-0.05/2)*0.8),
               interval_type = 'confidence',
               interval_level = 0.05);
    The false discovery risk is: 0.0594986

    False discovery risk with 1/8 likelihood support level

    >>> sgpv.risk(sgpval = 0,  null_lo = np.log(1/1.1),
                  null_hi = np.log(1.1), std_err = 0.8,
                  null_weights = 'Point',  null_space = 0,
               alt_weights = 'Uniform',
               alt_space = (2 + 1*norm.ppf(1-0.041/2)*0.8,
                            2 - 1*norm.ppf(1-0.041/2)*0.8),
               interval_type = 'likelihood',
               interval_level = 1/8);
    0.0505552

   With truncated normal weighting distribution

    >>> sgpv.risk(sgpval = 0,  null_lo = np.log(1/1.1),
                  null_hi = np.log(1.1), std_err = 0.8,
                  null_weights = 'Point',
               null_space = 0,  alt_weights = 'TruncNormal',
               alt_space = (2 + 1*norm.ppf(1-0.041/2)*0.8,
                            2 - 1*norm.ppf(1-0.041/2)*0.8),
               interval_type = 'likelihood',
               interval_level = 1/8);
    The false discovery risk is: 0.0490258

    False discovery risk with likelihood support intervall and wider null hypothesis

    >>> sgpv.risk(sgpval = 0,  null_lo = np.log(1/1.5),
                  null_hi = np.log(1.5), std_err = 0.8,
                  null_weights = 'Point',  null_space = 0,
               alt_weights = 'Uniform',
               alt_space = (2.5 + 1*norm.ppf(1-0.041/2)*0.8,
                            2.5 - 1*norm.ppf(1-0.041/2)*0.8),
               interval_type = 'likelihood',
               interval_level = 1/8);
    The false discovery risk is: 0.0168835

    False confirmation risk example

    >>> sgpv.risk(sgpval = 1,  null_lo = np.log(1/1.5),
                  null_hi = np.log(1.5), std_err = 0.15,
                  null_weights = 'Uniform',
               null_space = (0.01 + 1*norm.ppf(1-0.041/2)*0.15,
                             0.01 - 1*norm.ppf(1-0.041/2)*0.15),
               alt_weights = 'Uniform',
               alt_space = (np.log(1.5), 1.25*np.log(1.5)),
               interval_type = 'likelihood',
               interval_level = 1/8);
    The false confirmatory risk is: 0.0305952
    """
    import numpy as np
    from scipy import integrate
    from scipy.stats import norm
    from sgpv import sgpv

    # Convert inputs into arrays for easier handling
    null_space = np.asarray(null_space, dtype=np.float64)
    alt_space = np.asarray(alt_space, dtype=np.float64)

    # Warnings
    if sgpval not in [0, 1]:
        raise ValueError('sgpval must take a value of 0 or 1 to use sgpv.risk')

    if interval_type not in ['confidence', 'likelihood']:
        raise ValueError(
            "Parameter 'interval_type' must be one of the following: \n \
             'confidence' \n  'likelihood' \n  (credible/bayesian currently \
                                                not supported for sgpv.risk)")
    if not 0 < pi0 < 1:
        raise ValueError(
            'Values for option pi0 need to lie within the exclusive 0 - 1 interval.\
         A prior probability outside of this interval is not sensible.\
         The default value assumes that both hypotheses are equally likely.')

    if null_weights not in ['Point', 'Uniform', 'TruncNormal']:
        raise ValueError('Option null_weights must be one of the following:\
                         Point, Uniform or TruncNormal.')
    if alt_weights not in ['Point', 'Uniform', 'TruncNormal']:
        raise ValueError('Option alt_weights must be one of the following:\
                         Point, Uniform or TruncNormal.')
    # Relevant quantities
    # FDR = (1 + P(SGPV=0 | H1 ) / P(SGPV=0 | H0 ) *  P(H1) / P(H0) ) ^ (-1)
    fdr = None
    # FCR = (1 + P(SGPV=1 | H0 ) / P(SGPV=1 | H1 ) *  P(H0) / P(H1) ) ^ (-1)
    fcr = None

    p_sgpv_h0 = None   # `P_sgpv_H0` = P(SGPV=0 | H0 )
    p_sgpv_h1 = None   # `p_sgpv_h1` = P(SGPV=1 | H1 )

    if sgpval == 0:
        def power(x): return sgpv.power(
            true=x,
            null_lo=null_lo,
            null_hi=null_hi,
            std_err=std_err,
            interval_type=interval_type,
            interval_level=interval_level,
            no_print=True).poweralt
    if sgpval == 1:
        def power(x): return sgpv.power(
            true=x,
            null_lo=null_lo,
            null_hi=null_hi,
            std_err=std_err,
            interval_type=interval_type,
            interval_level=interval_level,
            no_print=True).powernull
    if null_lo == null_hi:
        if any(null_lo != null_space):
            print(
                'For a point indifference zone, specification of a different\
                null_space not permitted; null_space set to be ',
                round(null_lo, 2), '.')
        p_sgpv_h0 = power(x=null_lo)
    # interval null
    if null_lo != null_hi:
        # P_sgpv_H0 @ point (=type I error at null_space)
        if null_weights == "Point":
            if null_space.size != 1:
                raise ValueError('Null space must be a vector of len 1 when\
                        using a point null probability distribution.')
            p_sgpv_h0 = power(x=null_space)

        # p_sgpv_h0 averaged: check 'null_space' input
        if null_weights in ['Uniform', 'TruncNormal']:
            if null_space.size < 2:
                raise ValueError(
                    'Null space must not be a point to use averaging methods.')

            if null_space.size == 2:
                # truncate bounds to edge of null if null.spcae falls outside
                # indifference zone
                if (max(null_space) > null_hi) | (min(null_space) < null_lo):
                    print('Null space must be inside originally specified \
                    null hypothesis; at least one null space bound has been truncated.')
                    if max(null_space) > null_hi:
                        null_space[max(null_space)] = null_hi
                    if min(null_space) < null_lo:
                        null_space[min(null_space)] = null_lo

        # p_sgpv_h0 averaged uniformly
        if null_weights == 'Uniform':
            p_sgpv_h0 = 1 / (max(null_space) - min(null_space)) * integrate.quad(
                power, min(null_space), max(null_space))[0]

        # p_sgpv_h0 averaged using truncated normal as weighting distribution
        # function
        if null_weights == "TruncNormal":
            # default: mean of Normal distr at midpoint of null.space
            truncNorm_mu = np.mean(null_space)
            # default: std. dev of Normal distr same as assumed for
            # estimator
            truncNorm_sd = std_err

            def integrand(x): return power(x) * (norm.pdf(x,
                                                          truncNorm_mu,
                                                          truncNorm_sd) * (norm.cdf(max(null_space),
                                                                                    truncNorm_mu,
                                                                                    truncNorm_sd) - norm.cdf(min(null_space),
                                                                                                             truncNorm_mu,
                                                                                                             truncNorm_sd))**(- 1))

            p_sgpv_h0 = integrate.quad(
                integrand, min(null_space), max(null_space))[0]
    # calculate p_sgpv_h1
    # p_sgpv_h1 @ point
    if alt_weights == "Point":
        if alt_space.size != 1:
            raise ValueError(
                'alt space must be a vector of len 1 when using a point \
                    alternative probability distribution')
        if (alt_space >= null_lo) & (alt_space <= null_hi):
            raise ValueError(
                'alternative space must be outside of \
                the originally specified indifference zone')

        p_sgpv_h1 = power(x=alt_space)

    # p_sgpv_h1 averaged: check 'alt_space' input
    if alt_weights in ['Uniform', 'TruncNormal']:
        if alt_space.size < 2:
            raise ValueError(
                'alt space must not be a point to use averaging methods')

        if alt_space.size == 2:
            if np.all(alt_space > null_lo) & np.all(alt_space < null_hi):
                raise ValueError(
                    "Alternative space can not be contained inside \
                    indifference zone; 'null_space' and 'alt_space' \
                        might be flipped")
            if np.any(alt_space > null_lo) & np.any(alt_space < null_hi):
                raise ValueError(
                    "Alternative space can not intersect indifference zone")

    # p_sgpv_h1 averaged uniformly
    if alt_weights == 'Uniform':
        p_sgpv_h1 = 1 / (max(alt_space) - min(alt_space)) * integrate.quad(
            power, min(alt_space), max(alt_space))[0]

    # p_sgpv_h1 averaged using truncated normal as weighting
    # distribution function
    if alt_weights == "TruncNormal":
        # default: mean of Normal distr at midpoint of alt_space
        truncNorm_mu = np.mean(alt_space)
        # default: std. dev of Normal distr same as assumed for
        # estimator
        truncNorm_sd = std_err
        if np.any((truncNorm_mu, truncNorm_sd) is None):
            raise ValueError('trunNorm_mu` and `truncNorm_sd must be numeric;\
                             may not be None.')

        def integrand(x): return power(x) * (norm.pdf(x,
                                                      truncNorm_mu,
                                                      truncNorm_sd) * (norm.cdf(max(alt_space),
                                                                                truncNorm_mu,
                                                                                truncNorm_sd) - norm.cdf(min(alt_space),
                                                                                                         truncNorm_mu,
                                                                                                         truncNorm_sd)) ** (- 1))
        p_sgpv_h1 = integrate.quad(
            integrand, min(alt_space), max(alt_space))[0]

    # Calculate FDR or FCR
    if sgpval == 0:
        fdr = (1 + p_sgpv_h1 / p_sgpv_h0 * (1 - pi0) / pi0) ** (- 1)

    if sgpval == 1:
        fcr = (1 + p_sgpv_h0 / p_sgpv_h1 * pi0 / (1 - pi0)) ** (- 1)

    if fdr is not None:
        fdcr = round(fdr, 7)
        print('The false discovery risk (fdr) is:', fdcr)

    elif fcr is not None:
        fdcr = round(fcr, 7)
        print('The false confirmatory risk (fdr) is:', fdcr)

    return fdcr


def plot(*, est_lo, est_hi, null_lo, null_hi,
         set_order='sgpv', x_show=None, null_col=(0.815, 0.847, 0.909, 1),
         int_col=('cornflowerblue', 'darkslateblue', 'firebrick'),
         plot_axis=True, null_pt=None, outline_zone=True,
         title_lab="Title", x_lab="Position (by set_order)",
         y_lab="Outcome label", legend=True):
    """Plot interval estimates according to Second-Generation p-value rankings.

    This function displays user supplied interval estimates
    (support intervals, confidence intervals, credible intervals, etc.)
    according to its associated second-generation p-value ranking.

    Parameters
    ----------
    est_lo : array_like
        An array or list of lower bounds of interval estimates. Values must be
        finite for interval to be drawn. Must be of same length as est_hi.
    est_hi : array_like
        A array or list of upper bounds of interval estimates. Values must be
        finite for interval to be drawn. Must be of same length as est_lo.
    null_lo : float
        A scalar representing the lower bound of null interval (indifference zone).
        Value must be finite.
    null_hi : float
        A scalar representing the upper bound of null interval (indifference zone).
        Value must be finite.
    set_order : str or pandas.Series , optional
         A numeric vector/data series giving the desired order along the x-axis.
         If set_order is set to 'sgpv',
             the second-generation p-value ranking is used.
         If set_order is set to None,
             the original input ordering is used.
         If set_order is a pandas.
         Default is 'sgpv'.
    x_show : int, optional
         A scalar representing the maximum ranking on the x-axis that
         is displayed. Default is to display all intervals.
    null_col : valid_color, optional
        Coloring of the null interval (indifference zone).
        The default is Hawkes Blue: (0.815, 0.847, 0.909, 1) .
    int_col : valid_colors, optional
        Coloring of the intervals according to SGPV ranking.
        Default is ('cornflowerblue', 'darkslateblue', 'firebrick') for SGPVs
        of 1, between 0 and 1, and 0 respectively.
        Provide a list of three colors to replace the default colors.
    plot_axis : bool, optional
        Toggle for default axis plotting. The default is True.
    null_pt : float, optional
        A scalar representing a point null hypothesis. If set,
        the function will draw a horizontal dashed black line at this location.
        The default is None.
    outline_zone : bool, optional
       Toggle for drawing a slim white outline around the null zone.
       Helpful visual aid when plotting many intervals. The default is True.
    title_lab : str, optional
        Title text. The default is "Title".
    x_lab : str, optional
        x-axis label. The default is "Position (by setorder)".
    y_lab : str, optional
        y-axis label. The default is "Outcome label".
    legend : bool, optional
        Toggle for plotting the legend. The default is True.

    Return
    ------

    None.

    Example
    -------

    >>> from sgpv import sgpv
    >>> from sgpv import data
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> df = data.load()  # Load the example dataset as a dataframe
    >>> est_lo=df['ci.lo']
    >>> est_hi=df['ci.hi']
    >>> pvalue=df['p.value']
    >>> null_lo=-0.3
    >>> null_hi=0.3
    >>> title_lab="Leukemia Example"
    >>> y_lab="Fold Change (base 10)"
    >>> x_lab="Classical p-value ranking"
    >>> sgpv.plot(est_lo=est_lo, est_hi=est_hi, null_lo=null_lo,
    ...              null_hi=null_hi, set_order=pvalue, null_pt=0,
    ...              x_show=7000, outline_zone=True,
    ...            title_lab=title_lab, y_lab=y_lab, x_lab=x_lab );
    >>> plt.yticks(ticks=np.round(np.log10(np.asarray(
    ...        (1/1000,1/100,1/10,1/2,1,2,10,100,1000))),2), labels=(
    ...                           '1/1000','1/100','1/10',
                                    '1/2',1,2,10,100,1000));
    >>> plt.show()

    Note
    ----
    Some options and details of the R-code were not implemented because
    they do not exist in matplotlib or are difficult to create for only small gains.
   """

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib import collections as mc
    # import sgpv
    # Convert inputs
    # Create dataframe to make sorting and using x_show option easier
    data = pd.DataFrame({'est_lo': est_lo, 'est_hi': est_hi,
                         'null_lo': null_lo, 'null_hi': null_hi})
    est_lo = np.asarray(est_lo)
    est_hi = np.asarray(est_hi)

    # Errors -> Add checks for set_order, whether a variable in a dataset
    # exists or not
    if data['null_lo'].unique().size != data['null_hi'].unique().size:
        raise ValueError('null_lo and null_hi of different lengths')

    if data['est_lo'].unique().size != data['est_hi'].unique().size:
        raise ValueError('est_lo and est_hi of different lengths')

    if data['null_lo'].unique().size != 1 | data['null_hi'].unique().size != 1:
        raise ValueError('null_lo and null_hi must be scalars')

    if (not isinstance(set_order, (str, pd.core.series.Series))) and (
            set_order is not None):
        raise ValueError(
            'set_order can be only a string which equals "sgpv" or\
                          a pandas Series.')
    if isinstance(set_order, str) and set_order != 'sgpv':
        raise ValueError('set_order allows only "sgpv" as a string input.')

    # Set plot limits
    if x_show is None:
        x_show = len(est_lo)

    x_max = len(est_lo)
    x = np.arange(1, x_max + 1, 1)
    x_limits = (1, np.minimum(x_show, x_max))
    y_limits = (np.floor(np.minimum(np.amin(est_lo), np.amin(est_hi))),
                np.ceil(np.maximum(np.amax(est_lo), np.amax(est_hi))))
    # Compute SGPVs
    sgpvs = value(
        est_lo=est_lo,
        est_hi=est_hi,
        null_lo=null_lo,
        null_hi=null_hi)
    sgpvs.deltagap[np.where(sgpvs.pdelta == 0)] = - \
        sgpvs.deltagap[np.where(sgpvs.pdelta == 0)]
    data['pdelta'] = sgpvs.pdelta
    data['deltagap'] = sgpvs.deltagap
    sgpv_combo = np.where(sgpvs.pdelta == 0, sgpvs.deltagap, sgpvs.pdelta)
    sgpv_combo = np.array(sgpv_combo, dtype=np.float64)

    # Set order of x-axis
    if set_order is None:
        set_order = x
    if isinstance(set_order, str):
        if set_order == 'sgpv':
            set_order = sgpv_combo
    data['set_order'] = set_order

    # Sort arrays and then reducing size according to option x_show
    data = data.sort_values('set_order')
    # Subset intervals by SGPV value for coloring
    set_out = np.where(data['pdelta'] == 0)
    set_both = np.where((data['pdelta'] < 1) & (data['pdelta'] > 0))
    set_in = np.where(data['pdelta'] == 1)
    sets = [set_both, set_in, set_out]

    # Plotting
    fig, ax = plt.subplots()
    ax.set_xlabel(x_lab)
    ax.set_ylabel(y_lab)
    ax.set_title(title_lab)
    ax.set(xlim=x_limits, ylim=y_limits)

    # Not sure how to turn off only one axis
    if plot_axis is not True:
        ax.axis('off')

    labels = [r'$p_\delta$ = 0', r'0<$p_\delta$<1', r'$p_\delta$ = 1']
    # Null interval
    patches = []
    ax.fill_between(x, np.repeat(null_lo, len(x)), np.repeat(null_hi, len(x)),
                    color=null_col, label='Interval Null')
    patches.append(mpatches.Patch(color=null_col, label='Interval Null'))

    # SGPV-intervals -> Use lineCollection because matplotlib does not have
    # a ranged bar plot similar to Stata or a segment plot like R
    for i in range(len(sets)):
        intervals = _genpoints(x[sets[i]], data['est_lo'].iloc[sets[i]],
                               data['est_hi'].iloc[sets[i]])
        lc = mc.LineCollection(intervals, colors=int_col[i], linewidths=0.8)
        ax.add_collection(lc)
        patches.append(mpatches.Patch(color=int_col[i], label=labels[i]))

    # Detail indifference zone
    if null_pt is not None:
        y_0 = x * 0
        ax.plot(x, y_0, linestyle='--', color='k', linewidth=1)

    if outline_zone:
        ax.plot(x, np.repeat(null_lo, len(x)), x, np.repeat(null_hi, len(x)),
                color='w', linewidth=0.8, label='Interval Null')

    # Legend
    if legend:
        ax.legend(handles=patches, loc='upper right')

    return fig, ax


# Additonal helper function for plotting to shorten the code
def _genpoints(xlist, lbound, ubound):
    """Generate the points necessary to plot the lines showing the interval estimates.

    Parameters
    ----------
    xlist : array_like
        Array of values on the x-axis.
    lbound : array_like
        Lower bounds of the estimates.
    ubound : array_like
        Upper bounds of the estimates.

    Returns
    -------
    points : list
        List of start and end points for plotting the interval estimates.
    """
    xlist = xlist.tolist().copy()
    lbound = lbound.tolist().copy()
    ubound = ubound.tolist().copy()
    # Alternative approach -> shorter, seen in
    # https://pandas.pydata.org/docs/getting_started/10min.html Stack example
    xtemp = list(zip(*[xlist, xlist]))
    inttemp = list(zip(*[lbound, ubound]))
    points = []
    for xnew, ynew in zip(xtemp, inttemp):
        points.append(list(zip(xnew, ynew)))

    return points

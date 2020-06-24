"""False Discovery Risk for Second-Generation p-values."""


def fdrisk(*, sgpval: int = 0, null_lo: float, null_hi: float, std_err: float,
           interval_type: str, interval_level: float, null_weights: str, null_space,
           alt_weights: str, alt_space, pi0: float = 0.5):
    """False Discovery Risk for Second-Generation p-values.

    TODO: Suppress output from sgpower -> at the moment everything is printed
            in the console -> Python behaves different from R here
            Add user-defined expections for better error handling:
                see https://docs.python.org/3/tutorial/errors.html

    Parameters
    ----------
    sgpval : int, optional
        The observed second-generation p-value. Default is 0,
        which gives the false discovery risk.
    null_lo : float
        The lower bound of the indifference zone (null interval) upon which the
        second-generation p-value was based.
    null_hi : float
        The upper bound for the indifference zone (null interval) upon which the
        second-generation p-value was based.
    std_err : float
        Standard error of the point estimate.
    interval_type : STRING
        Class of interval estimate used. This determines the functional form of
        the power function. Options are 'confidence' for a (1-$\alpha$)100 %
        confidence interval and 'likelihood' for a 1/k likelihood support interval.
    interval_level : float
        Level of interval estimate. If inttype is 'confidence', the level is $\alpha$.
        If interval_type is 'likelihood', the level is 1/k (not k).
    null_weights : STRING
        Probability distribution for the null parameter space.
        Options are currently 'Point', 'Uniform', and 'TruncNormal'.
    null_space : array_like
        Support of the null probability distribution. If null_weights is 'Point',
        then null_space is a scalar. If null_weights is 'Uniform', then
        null_space is a vector of length two.
    alt_weights : str
        Probability distribution for the alternative parameter space.
        Options are currently 'Point', 'Uniform', and 'TruncNormal'.
        If null_weights is 'Uniform', then null_space is a vector of length two.
    alt_space : array_like
        Support of the alternative probability distribution. If null_weights is 'Point',
        then null_space is a scalar. If null_weights is 'Uniform', then
        null_space is a vector of length two.
    pi0 : float, optional
        Prior probability of the null hypothesis. The default is 0.5.
        This value can be only between 0 and 1 (exclusive).
        A prior probability outside of this interval is not sensible.
        The default value assumes that both hypotheses are equally likely.

    Raises
    ------
    ValueError
        Indicates that some value was outside of the expected range or
        outside of the accepted options.

    Returns
    -------
    res: float
        Numeric scalar representing the False discovery risk (FDR) or
        false confirmation risk (FCR) for the observed second-generation p-value.
        If sgpval = 0, the function returns false discovery risk (FDR).
        If sgpval = 1, the function returns false confirmation risk (FCR).

    Examples
    -------
    >>> from fdrisk import fdrisk
    >>> import numpy as np
    >>> from scipy.stats import norm
    # false discovery risk with 95% confidence level
    >>> fdrisk(sgpval = 0, null_lo = np.log(1/1.1), null_hi = np.log(1.1),
               std_err = 0.8, null_weights = 'Uniform',
               null_space = (np.log(1/1.1), np.log(1.1)), alt_weights = 'Uniform',
               alt_space = (2 + 1*norm.ppf(1-0.05/2)*0.8, 2 - 1*norm.ppf(1-0.05/2)*0.8),
               interval_type = 'confidence', interval_level = 0.05)
    0.0594986
    
    # false discovery risk with 1/8 likelihood support level
    >>> fdrisk(sgpval = 0,  null_lo = np.log(1/1.1), null_hi = np.log(1.1),
               std_err = 0.8,  null_weights = 'Point',  null_space = 0,
               alt_weights = 'Uniform',
               alt_space = (2 + 1*norm.ppf(1-0.041/2)*0.8, 2 - 1*norm.ppf(1-0.041/2)*0.8),
               interval_type = 'likelihood',  interval_level = 1/8)
    0.0505552

    ## with truncated normal weighting distribution
    >>> fdrisk(sgpval = 0,  null_lo = np.log(1/1.1), null_hi = np.log(1.1),
               std_err = 0.8,  null_weights = 'Point',
               null_space = 0,  alt_weights = 'TruncNormal',
               alt_space = (2 + 1*norm.ppf(1-0.041/2)*0.8, 2 - 1*norm.ppf(1-0.041/2)*0.8),
               interval_type = 'likelihood',  interval_level = 1/8)
    0.0490258
    
    # false discovery risk with LSI and wider null hypothesis
    >>> fdrisk(sgpval = 0,  null_lo = np.log(1/1.5), null_hi = np.log(1.5),
               std_err = 0.8,  null_weights = 'Point',  null_space = 0,
               alt_weights = 'Uniform',
               alt_space = (2.5 + 1*norm.ppf(1-0.041/2)*0.8, 2.5 - 1*norm.ppf(1-0.041/2)*0.8),
               interval_type = 'likelihood',  interval_level = 1/8)
    0.0168835
    
    # false confirmation risk example
    >>> fdrisk(sgpval = 1,  null_lo = np.log(1/1.5), null_hi = np.log(1.5),
               std_err = 0.15,  null_weights = 'Uniform',
               null_space = (0.01 + 1*norm.ppf(1-0.041/2)*0.15,
                             0.01 - 1*norm.ppf(1-0.041/2)*0.15),
               alt_weights = 'Uniform',
               alt_space = (np.log(1.5), 1.25*np.log(1.5)),
               interval_type = 'likelihood',  interval_level = 1/8)
    0.0305952

    """
    import numpy as np
    from scipy import integrate
    from scipy.stats import norm
    from sgpower import sgpower

    # Convert inputs into arrays for easier handling
    null_space = np.asarray(null_space, dtype=np.float64)
    alt_space = np.asarray(alt_space, dtype=np.float64)

    # Warnings
    if sgpval not in [0, 1]:
        raise ValueError('sgpval must take a value of 0 or 1 to use fdrisk')

    if interval_type not in ['confidence', 'likelihood']:
        raise ValueError(
            "Parameter 'interval_type' must be one of the following: \n \
             'confidence' \n  'likelihood' \n  (credible/bayesian currently \
                                                not supported for fdrisk)")
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
    Fdr = None
    # FCR = (1 + P(SGPV=1 | H0 ) / P(SGPV=1 | H1 ) *  P(H0) / P(H1) ) ^ (-1)
    Fcr = None

    P_sgpv_H0 = None   # `P.sgpv.H0` = P(SGPV=0 | H0 )
    P_sgpv_H1 = None   # `P.sgpv.H1` = P(SGPV=1 | H1 )

    if sgpval == 0:
        def power(x): return sgpower(true=x, null_lo=null_lo, null_hi=null_hi,
                                     std_err=std_err, interval_type=interval_type,
                                     interval_level=interval_level, 
                                     no_print=True).poweralt
    if sgpval == 1:
        def power(x): return sgpower(true=x, null_lo=null_lo, null_hi=null_hi,
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
        P_sgpv_H0 = power(x=null_lo)
    # interval null
    if null_lo != null_hi:
        # P.sgpv.H0 @ point (=type I error at null_space)
        if null_weights == "Point":
            if null_space.size != 1:
                raise ValueError('Null space must be a vector of len 1 when\
                        using a point null probability distribution.')
            P_sgpv_H0 = power(x=null_space)

        # P.sgpv.H0 averaged: check 'null_space' input
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

        # P.sgpv.H0 averaged uniformly
        if null_weights == 'Uniform':
            P_sgpv_H0 = 1 / (max(null_space) - min(null_space)) * integrate.quad(
                power, min(null_space), max(null_space))[0]

        # P.sgpv.H0 averaged using truncated normal as weighting distribution function
        if null_weights == "TruncNormal":
            # default: mean of Normal distr at midpoint of null.space
            truncNorm_mu = np.mean(null_space)
            # default: std. dev of Normal distr same as assumed for
            # estimator
            truncNorm_sd = std_err
            def integrand(x): return power(x) * (norm.pdf(x, truncNorm_mu,
                          truncNorm_sd) * (norm.cdf(max(null_space),
                            truncNorm_mu, truncNorm_sd) - norm.cdf(
                            min(null_space), truncNorm_mu, truncNorm_sd))**(- 1))

            P_sgpv_H0 = integrate.quad(integrand, min(null_space), max(null_space))[0]
            print(P_sgpv_H0)
    # calculate P.sgpv.H1
    # P.sgpv.H1 @ point
    if alt_weights == "Point":
        if alt_space.size != 1:
            raise ValueError(
                'alt space must be a vector of len 1 when using a point \
                    alternative probability distribution')
        if (alt_space >= null_lo) & (alt_space <= null_hi):
            raise ValueError(
                'alternative space must be outside of \
                the originally specified indifference zone')

        P_sgpv_H1 = power(x=alt_space)

    # P.sgpv.H1 averaged: check `alt.space` input
    if alt_weights in ['Uniform', 'TruncNormal']:
        if alt_space.size < 2:
            raise ValueError('alt space must not be a point to use averaging methods')

        if alt_space.size == 2:
            if np.all(alt_space > null_lo) & np.all(alt_space < null_hi):
                raise ValueError(
                    "Alternative space can not be contained inside \
                    indifference zone; 'null_space' and 'alt_space' \
                        might be flipped")
       # if np.any(alt_space > null_lo and alt_space < null_hi): The original approach
            if np.any(alt_space > null_lo) & np.any(alt_space < null_hi):
                raise ValueError("Alternative space can not intersect indifference zone")

    # P.sgpv.H1 averaged uniformly
    if alt_weights == 'Uniform':
        P_sgpv_H1 = 1 / (max(alt_space) - min(alt_space)) * integrate.quad(
            power, min(alt_space), max(alt_space))[0]
        # print('P_sgpv_H1 :',P_sgpv_H1)

    # P.sgpv.H1 averaged using truncated normal as weighting
    # distribution function
    if alt_weights == "TruncNormal":
        # default: mean of Normal distr at midpoint of alt.space
        truncNorm_mu = np.mean(alt_space)
        # default: std. dev of Normal distr same as assumed for
        # estimator
        truncNorm_sd = std_err
        if np.any((truncNorm_mu, truncNorm_sd) is None):
            raise ValueError('trunNorm_mu` and `truncNorm_sd must be numeric;\
                             may not be None.')

        def integrand(x): return power(x) * (norm.pdf(x, truncNorm_mu,
                              truncNorm_sd) * (norm.cdf(max(alt_space),
                            truncNorm_mu, truncNorm_sd) - norm.cdf(min(alt_space),
                            truncNorm_mu, truncNorm_sd)) ** (- 1))
        P_sgpv_H1 = integrate.quad(integrand, min(alt_space), max(alt_space))[0]
        # print(P_sgpv_H1)

    # Calculate FDR or FCR
    if sgpval == 0:
        Fdr = (1 + P_sgpv_H1 / P_sgpv_H0 * (1 - pi0) / pi0) ** (- 1)
        # print(Fdr)
    if sgpval == 1:
        Fcr = (1 + P_sgpv_H0 / P_sgpv_H1 * pi0 / (1 - pi0)) ** (- 1)
        # print(Fcr)

    # if FDR, Fcr is None, so c(Fdr, Fcr)=Fdr, and if FCR, Fdr is None,
    # risk = namedtuple('risk', 'Fdr, Fcr')
    # res = risk(Fdr, Fcr)
    # print(res)
    # return res
    if Fdr is not None:
        res = round(Fdr, 7)
        # val = 'Fdr'
    elif Fcr is not None:
        res = round(Fcr, 7)
        # val = 'Fcr :'
    return res

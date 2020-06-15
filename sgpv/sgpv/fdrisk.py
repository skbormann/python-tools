 
from pyensae.python2r_helper import make_tuple # Will be removed after I figured out the correct way
#Need to make available somehow the sgpower function
def fdrisk(*,sgpval=0, null_lo, null_hi, std_err,
                   interval_type, interval_level, pi0 = 0.5,
                   null_weights, null_space, alt_weights, alt_space):
    import numpy as np
    from scipy import integrate
    from scipy.stats import norm
    import sgpv
    from termcolor import colored
    """False Discovery Risk for Second-Generation p-values"""
       # Warnings
   # if ~  (set(sgpval) & set(make_tuple(0, 1))):
    #    sgpv.stop('sgpval must take a value of 0 or 1 to use fdrisk')

   # if ~  (set(interval_type) & set(make_tuple('confidence', 'likelihood'))):
   #     sgpv.stop("Parameter `interval.type` must be one of the following: \n  * confidence \n  * likelihood \n  (credible/bayesian not currently supported for fdrisk)")
   
           # Relevant quantities
    # FDR = (1 + P(SGPV=0 | H1 ) / P(SGPV=0 | H0 ) *  P(H1) / P(H0) ) ^ (-1)
    Fdr = None
    # FCR = (1 + P(SGPV=1 | H0 ) / P(SGPV=1 | H1 ) *  P(H0) / P(H1) ) ^ (-1)
    Fcr = None

    P_sgpv_H0 = None   # `P.sgpv.H0` = P(SGPV=0 | H0 )
    P_sgpv_H1 = None   # `P.sgpv.H1` = P(SGPV=1 | H1 )
    
    if sgpval == 0:
        power = lambda x: sgpv.sgpower(
                            true=x,
                            null_lo=null_lo,
                            null_hi=null_hi,
                            std_err=std_err,
                            interval_type=interval_type,
                            interval_level=interval_level)[0] # .power_alt
    if sgpval == 1:
        power = lambda x: sgpv.sgpower(true=x, null_lo=null_lo, null_hi=null_hi, 
                                  std_err=std_err, 
                                  interval_type=interval_type, 
                                  interval_level=interval_level)[2]
    if null_lo == null_hi:
        if any(null_lo != null_space):
            print(colored('for a point indifference zone, specification of a different `null.space` not permitted; `null.space` set to be ','red'),round(null_lo,2),'.')
        P_sgpv_H0 = power(x=null_lo)
            # interval null
    if null_lo != null_hi:

        # P.sgpv.H0 @ point (=type I error at null.space)
        if null_weights == "Point":

            if len(null_space) != 1:
                sgpv.sgpv.stop(
                    'null space must be a vector of len 1 when using a point null probability distribution')
            P_sgpv_H0 = power(x=null_space)

          # P.sgpv.H0 averaged: check `null.space` input
        if set(null_weights) & set(make_tuple(
                'Uniform', 'GBeta', 'TruncNormal')):

            if len(null_space) < 2:
                sgpv.stop('null space must not be a point to use averaging methods')
            if len(null_space) == 2:

                # truncate bounds to edge of null if null.spcae falls outside
                # indifference zone
                if max(null_space) > null_hi | min(null_space) < null_lo:

                    print(colored(
                        'null space must be inside originally specified null hypothesis; at least one null space bound has been truncated', 'red'))
                    if max(null_space) > null_hi:
                        null_space[max(null_space)] = null_hi
                    if min(null_space) < null_lo:
                        null_space[min(null_space)] = null_lo

              # P.sgpv.H0 averaged uniformly
            if null_weights == 'Uniform':

                P_sgpv_H0 = 1 / (max(null_space) - min(null_space)) * integrate.quad(
                    power, min(null_space), max(null_space))[0]

              # P.sgpv.H0 averaged using generalized beta as weighting
              # distribution function
            if null_weights == "GBeta":

                print(colored('placeholder for future implementation of Generalized Beta null probability distribution','red')
                    )

                P_sgpv_H0 = None

              # P.sgpv.H0 averaged using truncated normal as weighting
              # distribution function
            if null_weights == "TruncNormal":

                # default: mean of Normal distr at midpoint of null.space
                truncNorm_mu = np.mean(null_space)
                # default: std. dev of Normal distr same as assumed for
                # estimator
                truncNorm_sd = std_err

                integrand = lambda x:power(x) * (norm.pdf(x,
                                      truncNorm_mu,
                                      truncNorm_sd) * (norm.cdf(max(null_space),
                                                             truncNorm_mu,
                                                             truncNorm_sd) 
                                                            - norm.cdf(min(null_space),
                                                              truncNorm_mu,
                                                              truncNorm_sd)) ** (- 1))

                P_sgpv_H0 = integrate.quad(integrand, min(null_space),
                                           max(null_space))[0]

        # calculate P.sgpv.H1

          # P.sgpv.H1 @ point
        if alt_weights == "Point":

            if len(alt_space) != 1:
                sgpv.stop(
                    'alt space must be a vector of len 1 when using a point alternative probability distribution')
            if ((alt_space >= null_lo) & (alt_space <= null_hi)):
                sgpv.stop(
                    'alternative space must be outside of the originally specified indifference zone')

            P_sgpv_H1 = power(x=alt_space)

          # P.sgpv.H1 averaged: check `alt.space` input
        if set(alt_weights) & set(make_tuple(
                'Uniform', 'GBeta', 'TruncNormal')):

            if len(alt_space) < 2:
                sgpv.stop('alt space must not be a point to use averaging methods')
            if len(alt_space) == 2:

                if all(alt_space > null_lo) & all(alt_space < null_hi):
                    sgpv.stop(
                        "alternative space can not be contained inside indifference zone; `null.space` and `alt.space` might be flipped")
                if any(alt_space > null_lo & alt_space < null_hi):
                    sgpv.stop("alternative space can not intersect indifference zone")

              # P.sgpv.H1 averaged uniformly
            if alt_weights == 'Uniform':

                P_sgpv_H1 = 1 / (max(alt_space) - min(alt_space)) * integrate.quad(
                    power, min(alt_space), max(alt_space))[0] 

              # P.sgpv.H1 averaged using generalized beta as weighting
              # distribution function
            if alt_weights == "GBeta":

                warning(
                    'placeholder for future implementation of Generalized Beta null probability distribution')
                P_sgpv_H1 = None

              # P.sgpv.H1 averaged using truncated normal as weighting
              # distribution function
            if alt_weights == "TruncNormal":

                # default: mean of Normal distr at midpoint of alt.space
                truncNorm_mu = mean(alt_space)
                # default: std. dev of Normal distr same as assumed for
                # estimator
                truncNorm_sd = std_err

                if any(is_na(make_tuple(truncNorm_mu, truncNorm_sd))):
                    sgpv.stop(
                        'error: `trunNorm.mu` and `truncNorm.sd` must be NULL or numeric; may not be NA')

                integrand = lambda x:power(x) * (norm.pdf(x,
                      truncNorm_mu,
                      truncNorm_sd) * (norm.cdf(max(alt_space),
                                             truncNorm_mu,
                                             truncNorm_sd) 
                                            - norm.cdf(min(alt_space),
                                              truncNorm_mu,
                                              truncNorm_sd)) ** (- 1))
                P_sgpv_H1 = integrate.quad(
                    integrand,
                    min(alt_space),
                    max(alt_space))[0]

              # Calculate FDR or FCR
            if sgpval == 0:
                Fdr = (1 + P_sgpv_H1 / P_sgpv_H0 * (1 - pi0) / pi0) ** (- 1)
            if sgpval == 1:
                Fcr = (1 + P_sgpv_H0 / P_sgpv_H1 * pi0 / (1 - pi0)) ** (- 1)

            # if FDR, Fcr is null, so c(Fdr, Fcr)=Fdr, and if FCR, Fdr is null,
            # so c(Fdr, Fcr)=Fcr
            return (Fdr, Fcr)

        
        
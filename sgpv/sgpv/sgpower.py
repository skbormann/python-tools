# Add type hints for parameters
def sgpower(*,true, null_lo, null_hi, std_err=1, interval_type, interval_level, bonus=False):
    """ sgpower - Compute power/type I error for Second-Generation p-values approach """
    import numpy as np
    from scipy import integrate
    from scipy.stats import norm
    
    
    #Need to learn error handling
    #if !intervaltype in ['confidence', 'likelihood'] #stop("Parameter `intervaltype` must be one of the following: \n  * confidence \n  * likelihood \n  (credible not currently supported for power calculations)")

    if interval_type == 'confidence': 
        Z = norm.ppf(1-interval_level/2)
    if interval_type == 'likelihood':
        Z = norm.ppf(1-2*norm.cdf(-np.sqrt(2*np.log(1/interval_level)))/2)
  
    ## P(SGPV=0 | true )    (see Blume et al. (2018) eq.(S4) for CI/LSI)
    power0 = norm.cdf(null_lo/std_err - true/std_err - Z) + norm.cdf(-null_hi/std_err + true/std_err - Z)
  
    ## P(SGPV=1 | true )    (see Blume et al. (2018) eq.(S7) for CI/LSI)
    # -> only for symmetric null hypothesis
    if (null_hi-null_lo) >= 2*Z*std_err: 
        power1 = norm.cdf(null_hi/std_err - true/std_err - Z) - norm.cdf(null_lo/std_err - true/std_err + Z)
    
    if ((null_hi-null_lo) < 2*Z*std_err): 
        power1 = 0
    
  
    ## P(0<SGPV<1 | true)   (see Blume et al. (2018) eq.(S8, S9) for CI/LSI)
    # -> only for symmetric null hypothesis
    if ((null_hi-null_lo) <= 2*Z*std_err): 
        powerinc = 1 - norm.cdf(null_lo/std_err - true/std_err - Z) - norm.cdf(-null_hi/std_err + true/std_err - Z)
    
    if ((null_hi-null_lo) > 2*Z*std_err): 
        powerinc = 1 - (norm.cdf(null_lo/std_err - true/std_err - Z) + norm.cdf(-null_hi/std_err + true/std_err - Z)) - (norm.cdf(null_hi/std_err - true/std_err - Z) - norm.cdf(null_lo/std_err - true/std_err + Z))
    
  
  
    ## check
    #if (any(round(power0+powerinc+power1,7) != 1)) warning(paste0('error: power0+powerinc+power1 != 1 for indices ', paste(which(round(power0+powerinc+power1,7) != 1), collapse=', ')))
  
    ## bonus: type I error summaries
    pow0 = lambda x: norm.cdf(null_lo/std_err - x/std_err - Z) + norm.cdf(-null_hi/std_err + x/std_err - Z)
  
    minI = pow0((null_lo+null_hi)/2)
    maxI = pow0(null_lo)
    avgI = 1/(null_hi-null_lo)*integrate.quad(pow0, null_lo, null_hi)[0]
  
    typeI = (minI, maxI, avgI)
    if null_lo<=0 & 0<=null_hi:
       # print('at 0'=pow0(0), 'min'=minI, 'max'=maxI, 'mean'=avgI)
        #typeI = (pow0(0), minI, maxI, avgI)
        typeI = ('at 0 =', round(pow0(0),7), 'min =', round(minI,7), 'max =',round(maxI,7), 'mean =', round(avgI,7))
  
    ## more bonus: P(inc | null) (basically analogous to type II error but for confirmation)
    ## (TBD)
    #Displaying of bonus statistics not correct yet
   # print('poweralt' = power0, 'powerinc' = powerinc, 'powernull' = power1, 'type I error summaries' = typeI)
    print('poweralt =', round(power0,7), 'powerinc =', round(powerinc,7), 'powernull = ', round(power1,7), '\n type I error summaries: \n', typeI)
    #return('poweralt' : power0, 'powerinc' : powerinc, 'powernull' : power1, 'typeI' : typeI)
    return(power0, powerinc, power1, typeI)
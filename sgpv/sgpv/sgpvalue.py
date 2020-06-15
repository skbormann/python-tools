def sgpvalue(*, nulllo, nullhi, estlo, esthi, infcorrection=1e-5, warnings=True):
    """sgpvalue computes Second-Generation p-values and delta-gaps"""
    #Convert inputs into np.array to emulate R behaviour
    import numpy as np
    if len(nullhi) != len(nulllo):
        print("nulllo and nullhi are of different length")
    if (len(estlo) != len(esthi)):
        print('est.lo and est.hi of different lengths')

    #
    if len(nulllo) == 1:
        nulllo = np.repeat(nulllo, len(estlo))
        nullhi = np.repeat(nullhi, len(estlo))

    # Compute Interval Lengths
    estlen = np.array(esthi) - np.array(estlo)
    nulllen = np.array(nullhi) - np.array(nulllo)

    # Warnings -> to be added once I know how to check for these

    # SGPV computation
    overlap = np.minimum(esthi, nullhi) - np.maximum(estlo, nulllo)
    overlap = np.maximum(overlap, 0)

    bottom = np.minimum(2 * nulllen, estlen)

    pdelta = overlap / bottom

    # Zero-length & Infinite-length intervals -> to be added once I know how to check for these

    # Calculate delta gap
    # deltagap = np.repeat(NaN, len(pdelta))
    # deltagap[! is.na(p.delta) & (p.delta == 0)] = 0

    gap = np.maximum(estlo, nulllo) - np.minimum(nullhi, esthi)

    delta = nulllen / 2

    # Report unscaled delta gap if null has infinite length
    # delta[nulllen == Inf] = 1

    # Report unscaled delta gap if null has length zero
    # delta[nulllen == 0] = 1

    dg = gap / delta

    # deltagap[! is.na(pdelta) & (pdelta == 0)] = dg[! is.na(pdelta) & (pdelta == 0)]

    return pdelta  # , deltagap

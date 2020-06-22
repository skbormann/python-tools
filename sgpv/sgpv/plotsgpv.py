"""Plot interval estimates according to Second-Generation p-value rankings."""


def plotsgpv(*, est_lo, est_hi, null_lo, null_hi,
             set_order='sgpv', x_show=None, null_col=(0.815, 0.847, 0.909, 1),
             int_col=('cornflowerblue', 'darkslateblue', 'firebrick'),
             plot_axis=True, null_pt=None, outline_zone=True,
             title_lab="Title", x_lab="Position (by set_order)",
             y_lab="Outcome label", legend_on=True):
    """
    Plot interval estimates according to Second-Generation p-value rankings
    NOTE: Some options and details of the R-code were not implemented because
    they do not exist in matplotlib or are difficult to create for only small gains.

    Parameters
    ----------
    est_lo : array_like
        A numeric vector of lower bounds of interval estimates. Values must be
        finite for interval to be drawn. Must be of same length as 'est_hi'.
    est_hi : array_like
        A numeric vector of upper bounds of interval estimates. Values must be
        finite for interval to be drawn. Must be of same length as 'est_lo'.
    null_lo : float
        A scalar representing the lower bound of null interval (indifference zone).
        Value must be finite.
    null_hi : float
        A scalar representing the upper bound of null interval (indifference zone).
        Value must be finite.
    set_order : TYPE, optional
         The default is 'sgpv'. A numeric vector giving the desired order along the x-axis.
         If 'set_order' is set to 'sgpv', the second-generation p-value ranking is used.
         If 'set_order' is set to None, the original input ordering is used.
    x_show : int, optional
         A scalar representing the maximum ranking on the x-axis that
         is displayed. Default is to display all intervals.
    null_col : TYPE, optional
        Coloring of the null interval (indifference zone).
        The default is Hawkes Blue: (0.815, 0.847, 0.909, 1) .
    int_col : TYPE, optional
        Coloring of the intervals according to SGPV ranking.
        Default is ('cornflowerblue', 'darkslateblue', 'firebrick') for SGPVs
        of 1, in (0,1), and 0 respectively.
    plot_axis : TYPE, optional
        Toggle for default axis plotting. The default is True.
    null_pt : TYPE, optional
        A scalar representing a point null hypothesis. If set,
        the function will draw a horizontal dashed black line at this location.
        The default is None.
    outline_zone : TYPE, optional
        DESCRIPTION. The default is True.
    title_lab : TYPE, optional
        Title text. The default is "Title".
    x_lab : str, optional
        x-axis label. The default is "Position (by setorder)".
    y_lab : str, optional
        y-axis label. The default is "Outcome label".
    legend_on : str, optional
        Toggle for plotting the legend. The default is True.
    plot_opt : ,optional
        Additional plotting options. Only options for matplotlib are valid and
        will override previous settings.

    Returns
    -------
    None.


    Examples
    --------
    import pandas as pd
    df = pd.read_csv('../data/leukstats.csv', index_col=0)
    est_lo=df['ci.lo']
    est_hi=df['ci.hi']
    pvalue=df['p.value']
    null_lo=-0.3
    null_hi=0.3
    title_lab="Leukemia Example"
    y_lab="Fold Change (base 10)"
    x_lab="Classical p-value ranking"
    plotsgpv(est_lo=est_lo, est_hi=est_hi, null_lo=null_lo, null_hi=null_hi,
             set_order=pvalue, null_pt=0, x_show=7000, outline_zone=True,
             title_lab=title_lab, y_lab=y_lab, x_lab=x_lab )
    plt.yticks(ticks=np.round(np.log10(np.asarray(
        (1/1000,1/100,1/10,1/2,1,2,10,100,1000))),2), labels=(
                           '1/1000','1/100','1/10','1/2',1,2,10,100,1000))
   plt.show()
   """

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from sgpvalue import sgpvalue

    # Convert inputs
    # Create dataframe to make sorting and using x_show option easier
    data = pd.DataFrame({'est_lo': est_lo, 'est_hi': est_hi, 'null_lo': null_lo,
                         'null_hi': null_hi})
    est_lo = np.asarray(est_lo)
    est_hi = np.asarray(est_hi)

    # Errors

    if data['null_lo'].unique().size != data['null_hi'].unique().size:
        raise ValueError('null_lo and null_hi of different lengths')

    if data['est_lo'].unique().size != data['est_hi'].unique().size:
        raise ValueError('est_lo and est_hi of different lengths')

    if data['null_lo'].unique().size != 1 | data['null_hi'].unique().size != 1:
        raise ValueError('null_lo and null_hi must be scalars')

    # Set plot limits
    if x_show is None:
        x_show = len(est_lo)

    x_max = len(est_lo)
    x = np.arange(1, x_max+1, 1)  # Is x really needed or just a take over from R?
    x_limits = (1, np.minimum(x_show, x_max))  
    y_limits = (np.floor(np.minimum(np.amin(est_lo), np.amin(est_hi))),
                np.ceil(np.maximum(np.amax(est_lo), np.amax(est_hi))))

    # Compute SGPVs
    sgpv = sgpvalue(est_lo=est_lo, est_hi=est_hi, null_lo=null_lo, null_hi=null_hi)
    sgpv.deltagap[np.where(sgpv.pdelta == 0)] = -sgpv.deltagap[np.where(sgpv.pdelta == 0)]
    data['pdelta'] = sgpv.pdelta
    data['deltagap'] = sgpv.deltagap
    sgpv_combo = np.where(sgpv.pdelta == 0, sgpv.deltagap, sgpv.pdelta)
    sgpv_combo = np.array(sgpv_combo, dtype=np.float64)
    # data['sgpv_combo'] = sgpv_combo

    # Set order of x-axis
    if set_order is None:
        set_order = x
    if set_order is 'sgpv':
        set_order = sgpv_combo
    data['set_order'] = set_order

    # Sort arrays and then reducing size according to option x_show
    #  pdelta = sgpv.pdelta[set_order].copy()
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

    labels = ['$p_\delta$ = 0', '0<$p_\delta$<1', '$p_\delta$ = 1']
    # Null interval
    patches = []
    ax.fill_between(x, np.repeat(null_lo, len(x)), np.repeat(null_hi, len(x)),
                    color=null_col, label='Interval Null')
    patches.append(mpatches.Patch(color=null_col, label='Interval Null'))

    # SGPV-intervals
    # Color is not correctly set
    for i in range(len(sets)):
        interval = genpoints(x[sets[i]], data['est_lo'].iloc[sets[i]],
                             data['est_hi'].iloc[sets[i]])
        intervals(interval, int_col[i], ax)
        patches.append(mpatches.Patch(color=int_col[i], label=labels[i]))

    # Detail indifference zone
    if null_pt is not None:
        y_0 = x*0
        ax.plot(x, y_0, linestyle='--', color='k', linewidth=1)

    if outline_zone:
        ax.plot(x, np.repeat(null_lo, len(x)), x, np.repeat(null_hi, len(x)),
                color='w', linewidth=0.8, label='Interval Null')

    # Legend
    if legend_on:
        ax.legend(handles=patches, loc='upper right')

    return fig, ax


def genpoints(xlist, lbound, ubound):
    """
    Parameters
    ----------
    xlist : TYPE
        DESCRIPTION.
    lbound : TYPE
        DESCRIPTION.
    ubound : TYPE
        DESCRIPTION.

    Returns
    -------
    points : TYPE
        DESCRIPTION.

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


def intervals(points, color, figure):
    """

    Parameters
    ----------
    points : TYPE
        DESCRIPTION.
    color : TYPE
        DESCRIPTION.
    figure : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    from matplotlib import collections as mc

    lc = mc.LineCollection(points, colors=color, linewidths=0.8)
    figure.add_collection(lc)

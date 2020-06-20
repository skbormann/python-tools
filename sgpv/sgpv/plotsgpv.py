def plotsgpv(*, est_lo, est_hi, null_lo, null_hi,
             set_order='sgpv', x_show=None, null_col=(0.815, 0.847, 0.909, 1),
             int_col=('firebrick', 'cornflowerblue', 'darkslateblue'),
             plot_axis=True, null_pt=None, outline_zone=True,
             title_lab="Title", x_lab="Position (by setorder)", y_lab="Outcome label",
             legend_on=True):
    """
    Plot ...
    To-Do: Sort does not work and color for pdelta=0 does not work

    Parameters
    ----------
    * : TYPE
        DESCRIPTION.
    est_lo : TYPE
        DESCRIPTION.
    est_hi : TYPE
        DESCRIPTION.
    null_lo : TYPE
        DESCRIPTION.
    null_hi : TYPE
        DESCRIPTION.
    set_order : TYPE, optional
        DESCRIPTION. The default is 'sgpv'.
    x_show : TYPE, optional
        DESCRIPTION. The default is NaN.
    #null_col : TYPE, optional
        DESCRIPTION. The default is rgb(208,216,232,max=255).
    #int_col : TYPE, optional
        DESCRIPTION. The default is c("cornflowerblue","firebrick3","darkslateblue").
    plot_axis : TYPE, optional
        DESCRIPTION. The default is [True,True].
    null_pt : TYPE, optional
        DESCRIPTION. The default is NA.
    outline_zone : TYPE, optional
        DESCRIPTION. The default is True.
    title_lab : TYPE, optional
        DESCRIPTION. The default is "Title".
    x_lab : str, optional
        DESCRIPTION. The default is "Position (by setorder)".
    y_lab : str, optional
        DESCRIPTION. The default is "Outcome label".
    legend_on : str, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    None.


    Examples
    --------
    import pandas as pd
    df = pd.read_csv('../data/leukstats.csv', index_col=1)
    est_lo=df['ci.lo'].copy()
    est_hi=df['ci.hi'].copy()
    pvalue=df['p.value'].copy()
    null_lo=-0.3
    null_hi=0.3
    title_lab="Leukemia Example"
    y_lab="Fold Change (base 10)"
    x_lab="Classical p-value ranking"
    plotsgpv(est_lo=est_lo, est_hi=est_hi, null_lo=null_lo, null_hi=null_hi,
             set_order=pvalue.to_numpy(), title_lab=title_lab, y_lab=y_lab, x_lab=x_lab)

    """

    import matplotlib.pyplot as plt
    #from matplotlib import collections as mc
    import matplotlib.patches as mpatches
    import numpy as np
    from sgpvalue import sgpvalue
    from stop import stop

    #Convert inputs
    est_lo = np.asarray(est_lo).copy()
    est_hi = np.asarray(est_hi).copy()

    #### Errors
    if len([null_lo]) != len([null_hi]):
        stop('null_lo and null_hi of different lengths')

    if len(est_lo) != len(est_hi):
        stop('est_lo and est_hi of different lengths')

    if len([null_lo]) != 1 | len([null_lo]) != 1:
        stop('null_lo and null_hi must be scalars')

    #### Set plot limits
    if x_show is None:
        x_show = len(est_lo)

    x_max = len(est_lo)
    x = np.arange(1, np.minimum(x_show, x_max)+1, 1)

    x_limits = (1, np.minimum(x_show, x_max)) # no removal of None yet
    y_limits = (np.floor(np.minimum(np.amin(est_lo), np.amin(est_hi))), \
             np.ceil(np.maximum(np.amax(est_lo), np.amax(est_hi))))

    #### Compute SGPVs
    sgpv = sgpvalue(est_lo=est_lo, est_hi=est_hi, null_lo=null_lo, null_hi=null_hi)
    sgpv.deltagap[np.where(sgpv.pdelta == 0)] = -sgpv.deltagap[np.where(sgpv.pdelta == 0)]
    sgpv_combo = np.where(sgpv.pdelta == 0, sgpv.deltagap, sgpv.pdelta)
    sgpv_combo = np.array(sgpv_combo, dtype=np.float64)

    #### Set order of x-axis -> not correct yet -> taking care of option x_show not working yet
    if set_order is None:
        set_order = x.argsort()
    elif set_order == 'sgpv':
        set_order = np.argsort(sgpv_combo)
    elif set_order is not None:
        set_order = np.argsort(set_order)

    #Sort arrays and then reducing size according to option x_show 
    pdelta = sgpv.pdelta[set_order].copy()
    #### Subset intervals by SGPV value for coloring
#     	gap.marker <- 1*is.na(sgpv$delta.gap[set.order])

# 	set.out  <- x[gap.marker[x]==0]
# 	set.in   <- x[(gap.marker[x]==1) & (sgpv$p.delta[set.order]==1)]
# 	set.both <- x[(gap.marker[x]==1) & (sgpv$p.delta[set.order]<1)]
    
    #gap_marker = np.where(sgpv.deltagap[set_order] == None, 1,0)
# Not correctly setting the case of 0< pdelta <1
    set_out = np.where(pdelta[:len(x)] == 0)
    set_both = np.where((pdelta[:len(x)] < 1) & (pdelta[:len(x)]>0))
    set_in = np.where(pdelta[:len(x)] == 1)
    sets = [set_out, set_both, set_in]

    #### Plotting
    fig, ax = plt.subplots()
    ax.set_xlabel(x_lab)
    ax.set_ylabel(y_lab)
    ax.set_title(title_lab)
    ax.set(xlim=x_limits, ylim=y_limits)

    #Not sure how to turn off only one axis
    if plot_axis is not True:
        ax.axis('off')

    labels = ['$p_\delta$ = 0', '0<$p_\delta$<1', '$p_\delta$ = 1']
    #Null interval
    patches = []
    ax.fill_between(x, np.repeat(null_lo, len(x)), np.repeat(null_hi, len(x)),
                    color=null_col, label='Interval Null')
    patches.append(mpatches.Patch(color=null_col, label='Interval Null'))

    #This code contains a bug which prevents the correct calculation of the intervals for pdelta=0
    #SGPV-intervals
    #Color is not correctly set
    for i in range(len(sets)):
        interval = genpoints(x[sets[i]], est_lo[sets[i]], est_hi[sets[i]])
        intervals(interval, int_col[i], ax)
        patches.append(mpatches.Patch(color=int_col[i], label=labels[i]))

    #Detail indifference zone
    if null_pt is not None:
        y_0 = x*0
        ax.plot(x, y_0, linestyle='--', color='k', linewidth=1)

    if outline_zone:
        ax.plot(x, np.repeat(null_lo, len(x)), x, np.repeat(null_hi, len(x)),
                color='w', linewidth=0.8, label='Interval Null') # Need to change linewidth

    #Legend
    if legend_on:
        ax.legend(handles=patches, loc='upper right')


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

    # xtemp = []
    # for i in range(len(xlist)):
    #     xtemp.append((xlist[i], xlist[i]))
    
    #Alternative approach -> shorter, seen in https://pandas.pydata.org/docs/getting_started/10min.html Stack example
    xtemp = list(zip(*[xlist, xlist]))

    # inttemp = []
    # for i in range(len(xlist)):
    #     inttemp.append((lbound[i], ubound[i]))
    
    #Alternative approach
    inttemp = list(zip(*[lbound, ubound]))

    points = []
    for xnew, ynew in zip(xtemp, inttemp): #-> can be made shorter
        points.append(list(zip(xnew, ynew)))
    #Alternative approach
    #points = list(zip(*[xtemp, inttemp])) 
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
    from matplotlib import collections  as mc

    lc = mc.LineCollection(points, colors=color, linewidths=0.8)
    figure.add_collection(lc)

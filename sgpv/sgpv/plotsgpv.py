def plotsgpv(*,est_lo, est_hi, null_lo, null_hi,
        set_order="sgpv", x_show=NaN, #null_col=rgb(208,216,232,max=255),
        #int_col=c("cornflowerblue","firebrick3","darkslateblue"),
          plot_axis=[True,True],
        null_pt=None, outline_zone=True,
        title_lab="Title", x_lab="Position (by setorder)", y_lab="Outcome label",
        legend_on=True ):
    """
    

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
        DESCRIPTION. The default is "sgpv".
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
    df = pd.read_csv('../data/leukstats.csv', index_col=0)
    est_lo=df['ci.lo']
    est_hi=df['ci.hi']
    null_lo=-0.3
    null_hi=0.3
    title_lab="Leukemia Example"
    y_lab="Fold Change (base 10)"
    x_lab="Classical p-value ranking"

    """

    import matplotlib.pyplot as plt
    from matplotlib import collections  as mc
    import matplotlib.patches as mpatches
    import numpy as np
    from sgpvalue import sgpvalue
    from stop import stop
    
    #Color settings -> translated R-colors into RGB for Python \
        #-> Not sure how to install the colours in Python for easier referencing.

    # firebrick3 = (205, 38, 38)
    # cornflowerblue = (100, 149, 237)
    # darkslateblue = (72, 61, 139)


    intcoldefault = ('firebrick','cornflowerblue', 'darkslateblue' ) 
    nullcoldefault = (0.815, 0.847, 0.909,1) #Hawkes Blue
    #Convert inputs
    est_lo = np.asarray(est_lo)
    est_hi = np.asarray(est_hi)
    
    #### Errors
    if (len([null_lo])!=len([null_hi])) :
      stop('null_lo and null_hi of different lengths')
    

    if(len(est_lo)!=len(est_hi)) :
      stop('est_lo and est_hi of different lengths')
    

    if(len([null_lo])!=1 | len([null_lo])!=1) :
      stop('null_lo and null_hi must be scalars')
    

    #### Set plot limits
    x_max     = len(est_lo)
    x         = np.arange(1,np.minimum(x_show,x_max)+1,1)

    x_limits = (1,np.minimum(x_show,x_max)) # no removal of None yet
    y_limits = (np.floor(np.minimum(np.amin(est_lo),np.amin(est_hi))), \
             np.ceil(np.maximum(np.amax(est_lo),np.amax(est_hi))))

    #### Compute SGPVs
    sgpv = sgpvalue(est_lo=est_lo, est_hi=est_hi, null_lo=null_lo, null_hi=null_hi)
    sgpv_combo = np.where(sgpv.pdelta==0,-(sgpv.deltagap),sgpv.pdelta) #Not correct dealings yet wit None
    
    #### Set order of x-axis
    if set_order[1] == None:
        set_order = x
    if (set_order[1]=="sgpv"):
        set_order = np.sort(sgpv_combo)

    #### Subset intervals by SGPV value for coloring
    # gap_marker = 1*is.na(sgpv.deltagap[set_order])

    # set_out  = x[gap_marker[x]==0]
    # set_in   = x[(gap_marker[x]==1) & (sgpv.pdelta[set_order]==1)]
    # set_both = x[(gap_marker[x]==1) & (sgpv.pdelta[set_order]<1)]
    
    set_out = np.where(sgpv.pdelta[:len(x)]==0)  
    set_both = np.where(sgpv.pdelta[:len(x)]<1)
    set_in = np.where(sgpv.pdelta[:len(x)]==1)
    sets=[set_out, set_both, set_in]
    
    
    
    fig, ax =plt.subplots()
    ax.set_xlabel(x_lab)
    ax.set_ylabel(y_lab)
    ax.set_title(title_lab)
    ax.set(xlim=x_limits, ylim=y_limits)
    
#         #### Plotting
#     plot(1, null.pt, ylim=y.limits, xlim=x.limits, type="n", yaxt="n", xaxt="n",
#         ylab=y.lab, xlab=x.lab, main=title.lab)

#     if (plot_axis[1]=="True"): #Not sure how to turn off only one axis
    #{axis(side=1)}
#     if (plot_axis[2]=="True"):
    #{axis(side=2,las=2)}
    #ax.axis('on')

#     rect(1, null.lo, x.max, null.hi, col=null.col, border=NA)
    labels =['$p_\delta$ = 0', '0<$p_\delta$<1', '$p_\delta$ = 1']
    #Null interval
    patches=[]
    ax.fill_between(x,np.repeat(null_lo,len(x)), np.repeat(null_hi,len(x)) \
                            ,color=nullcoldefault ,label='Interval Null')      
    patches.append(mpatches.Patch(color=nullcoldefault, label='Interval Null'))  
    
    for set in range(len(sets)):
        points = genpoints(x[sets[set]], est_lo[sets[set]], est_hi[sets[set]])
        intervals(points, intcoldefault[set])
        patches.append(mpatches.Patch(color=intcoldefault[set], label=labels[set]))
        
#     ## Intervals where 0<SGPV<1
     
#     points(x[set.both],est.lo[set.order][set.both],cex=int.cex,pch=int.pch,col=int.col[1])
#     points(x[set.both],est.hi[set.order][set.both],cex=int.cex,pch=int.pch,col=int.col[1])
#     segments(x[set.both], est.lo[set.order][set.both], x[set.both], est.hi[set.order][set.both], lty=1, col=int.col[1])

#     ## Intervals where SGPV==1
#     points(x[set.in],est.lo[set.order][set.in],cex=int.cex,pch=int.pch,col=int.col[3])
#     points(x[set.in],est.hi[set.order][set.in],cex=int.cex,pch=int.pch,col=int.col[3])
#     segments(x[set.in], est.lo[set.order][set.in], x[set.in], est.hi[set.order][set.in], lty=1, col=int.col[3])

#     ## Intervals where SGPV==0
      #generate the points for the linecollection
          
#     points(x[set.out],est.lo[set.order][set.out],cex=int.cex,pch=int.pch,col=int.col[2])
#     points(x[set.out],est.hi[set.order][set.out],cex=int.cex,pch=int.pch,col=int.col[2])
#     segments(x[set.out], est.lo[set.order][set.out], x[set.out], est.hi[set.order][set.out], lty=1, col=int.col[2])
    
#     ## Detail indifference zone
    if null_pt is not None:
        y0 = x*0
        ax.plot(x,y0, linestyle='--', color='k', linewidth=1) # Need to add linestyle and color and linewidth
#     abline(h=null.pt,lty=2)
    if outline_zone==True:
        ax.plot(x,np.repeat(null_lo,len(x)),x ,np.repeat(null_hi,len(x)),color='w', linewidth=0.8, label='Interval Null') # Need to change linewidth
           
#     #### Legend        
    if (legend_on==True):
        ax.legend(handles=patches,loc='upper right')
#     legend("topright",c("Interval Null", expression("p"[delta]*" = 0"),
#                         expression("0 < p"[delta]*" < 1"), expression("p"[delta]*" = 1")),
#             lty=1,col=c(null.col,int.col[2],int.col[1],int.col[3]),lwd=c(6,1.5,1.5,1.5),bty="n")
#     plt.show()                    
    
    
def genpoints(xlist, lbound,ubound):
    xlist=xlist.tolist()
    lbound=lbound.tolist()
    ubound=ubound.tolist()
    
    xtemp=[]
    for i in range(len(xlist)):
        xtemp.append((xlist[i],xlist[i]))
    
    inttemp=[]
    for i in range(len(xlist)):
        inttemp.append((lbound[i],ubound[i]))
   
    points=[]
    for x,y in zip(xtemp,inttemp):
        points.append(list(zip(x,y)))
    
    return (points)

def intervals(points,color):
    lc = mc.LineCollection(points, colors=color,linewidths=0.8)
    ax.add_collection(lc)
    
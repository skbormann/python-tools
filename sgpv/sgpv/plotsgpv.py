def plotsgpv(*,est_lo, est_hi, null_lo, null_hi,
        set_order="sgpv", x_show=NaN, #null_col=rgb(208,216,232,max=255),
        #int_col=c("cornflowerblue","firebrick3","darkslateblue"),
          plot_axis=[True,True],
        null_pt=NA, outline_zone=True,
        title_lab="Title", x_lab="Position (by setorder)", y_lab="Outcome label",
        legend_on=True ):
    """plotsgpv """
    import matplotlib as plt
    import numpy as np
    import sgpvalue, stop
    
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
    x         = np.arange(1,x_max+1,1)

    x_limits = (1,np.minimum(x_show,x_max)) # no removal of None yet
    y_limits = (np.floor(np.minimum(np.amin(est_lo),np.amin(est_hi))), \
             np.ceil(np.maximum(np.amax(est_lo),np.amax(est_hi))))

    #### Compute SGPVs
    sgpv = sgpvalue(est_lo=est_lo, est_hi=est_hi, null_lo=null_lo, null_hi=null_hi)
    sgpv_combo = np.where(sgpv['pdelta']==0,-(sgpv['deltagap']),sgpv['pdelta']) #Not correct dealings yet wit None
    
    #### Set order of x-axis
    if set_order[1] == None:
        set_order = x
    if (set_order[1]=="sgpv"):
        set_order = order(sgpv_combo)

    #### Subset intervals by SGPV value for coloring
    gap_marker = 1*is.na(sgpv$delta.gap[set.order])

    set_out  = x[gap.marker[x]==0]
    set_in   = x[(gap.marker[x]==1) & (sgpv$p.delta[set.order]==1)]
    set_both = x[(gap.marker[x]==1) & (sgpv$p.delta[set.order]<1)]
    
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

#     ## Intervals where 0<SGPV<1
#     points(x[set.both],est.lo[set.order][set.both],cex=int.cex,pch=int.pch,col=int.col[1])
#     points(x[set.both],est.hi[set.order][set.both],cex=int.cex,pch=int.pch,col=int.col[1])
#     segments(x[set.both], est.lo[set.order][set.both], x[set.both], est.hi[set.order][set.both], lty=1, col=int.col[1])

#     ## Intervals where SGPV==1
#     points(x[set.in],est.lo[set.order][set.in],cex=int.cex,pch=int.pch,col=int.col[3])
#     points(x[set.in],est.hi[set.order][set.in],cex=int.cex,pch=int.pch,col=int.col[3])
#     segments(x[set.in], est.lo[set.order][set.in], x[set.in], est.hi[set.order][set.in], lty=1, col=int.col[3])

#     ## Intervals where SGPV==0
#     points(x[set.out],est.lo[set.order][set.out],cex=int.cex,pch=int.pch,col=int.col[2])
#     points(x[set.out],est.hi[set.order][set.out],cex=int.cex,pch=int.pch,col=int.col[2])
#     segments(x[set.out], est.lo[set.order][set.out], x[set.out], est.hi[set.order][set.out], lty=1, col=int.col[2])

#     ## Detail indifference zone
#     abline(h=null.pt,lty=2)

#     if (outline.zone==TRUE) {
#             abline(h=null.lo, col="white", lty=1, lwd=0.8)
#             abline(h=null.hi, col="white", lty=1, lwd=0.8)
#                             }
#     #### Legend
#     if (legend_on==True):
    # ax.legend((line1,line2,line3,line4),('Interval Null', 'pdelta = 0', '0<pdelta<1', 'pdelta = 1'),loc='upper right')
#     legend("topright",c("Interval Null", expression("p"[delta]*" = 0"),
#                         expression("0 < p"[delta]*" < 1"), expression("p"[delta]*" = 1")),
#             lty=1,col=c(null.col,int.col[2],int.col[1],int.col[3]),lwd=c(6,1.5,1.5,1.5),bty="n")
#     plt.show()                    
    
    
    
    
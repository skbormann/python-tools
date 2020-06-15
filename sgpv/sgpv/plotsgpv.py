def plotsgpv(*,estlo, esthi, nulllo, nullhi,
		setorder="sgpv", xshow=NA, #nullcol=rgb(208,216,232,max=255),
		#intcol=c("cornflowerblue","firebrick3","darkslateblue"),
		  plotaxis=[True,True],
		nullpt=NA, outlinezone=True,
		titlelab="Title", xlab="Position (by setorder)", ylab="Outcome label",
		legendon=True ):
    """plotsgpv """
    import matplotlib as plt
    
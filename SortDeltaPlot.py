
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from sys import argv
import os, glob, time
from SortTimePlot import get_pars, get_headinfo,\
	flatline, SSS_theo, SSW_theo, SSt_theo
from ErrorRate import errorplot

##=============================================================================
def main():
	"""
	NAME
		SortDeltaPlot.py
		
	PURPOSE
	
	EXECUTION
		python SortDeltaPlot.py [dirname]
		
	OUTPUT
		A single png with four subplots, saved in data directory
	
	EXAMPLE
		python SortDeltaPlot.py Rsults/Sorting/
		
	NOTES
		- There must be an equal number of Hopfield and notfield files in
		the directory.
		- Assumed that all filesshare the same k-values (othr than Deltas and
		zeroes).
	
	HISTORY
		2016/01/28	Started
		2016/02/08	Notfield added
	"""
	me = "SortDeltaPlot.main: "
	t0 = time.time()
	
	verbose = True
	
	try: argv[1]
	except IndexError:
		print main.__doc__
		raise IndexError(me+"Please read docstring")
	
	## Outfile name and number of points to plot
	plotfile = argv[1]+"/DeltaPlots.png"
	npts = 500

	##-------------------------------------------------------------------------

	filelist = np.sort(glob.glob(argv[1]+"/*.txt"))
	numfiles = len(filelist)
	if numfiles%2 != 0:
		raise IOError(me+"Expecting an even number of files (Hopfield and Notfield).")
	## Separate by Hop and Not
	hoplist = filelist[:numfiles/2]
	notlist = filelist[numfiles/2:]
	
	Delta = np.zeros(numfiles)
	S_fin = np.zeros(numfiles)
	t_SS = np.zeros(numfiles)
	t_SS_th = np.zeros(numfiles)
	W_srt = np.zeros(numfiles)
	Wdot_SS = np.zeros(numfiles)
	
	## Get data from all files
	for i in range(numfiles):
	
		Delta[i] = get_pars(filelist[i])
		k = get_headinfo(filelist[i])
	
		data = np.loadtxt(filelist[i], skiprows=10, unpack=True)
		data = data[:, ::int(data.shape[1]/npts)]
		
		t, ent, work, trans_C, trans_I = data[[0,5,6,8,9]]
		N = int(data[[1,2,3,4],0].sum())
		del data
		
		ent /= N*np.log(2)	
		Dent, SSidx = flatline(ent, N, Delta, kerfrac=100)
		
		## Assuming entropy is flat and work is linear
		
		S_fin[i] = np.mean(ent[SSidx:])
		t_SS[i] = t[SSidx]
		t_SS_th[i] = SSt_theo(k)
		# Wdot_srt[i] = np.mean(work[:SSidx-int(npts/20)])/t[SSidx]
		W_srt[i] = work[SSidx]
		Wdot_SS[i] = np.mean(work[SSidx:]-work[SSidx])/(t[-1]-t[SSidx])
	
	## ----------------------------------------------------
	
	## Separate into H and N
	newshape = (2,numfiles/2)
	Delta = Delta.reshape(newshape)
	S_fin = S_fin.reshape(newshape)
	t_SS = t_SS.reshape(newshape)
	t_SS_th = t_SS_th.reshape(newshape)
	W_srt = W_srt.reshape(newshape)
	Wdot_SS  = Wdot_SS.reshape(newshape)
	
	## Sort by Delta
	sortind = np.argsort(Delta)
	Delta = Delta[:,sortind[0]]
	try:
		assert np.all(Delta[0]==Delta[1])
	except AssertionError:
		raise IOError(me+"Check files.\n"+filelist.tostring())
	S_fin = S_fin[:,sortind[0]]
	t_SS = t_SS[:,sortind[0]]
	t_SS_th = t_SS_th[:,sortind[0]]
	W_srt = W_srt[:,sortind[0]]
	Wdot_SS = Wdot_SS[:,sortind[0]]

	## ----------------------------------------------------
	
	## Get k values, assuming the same for all files in directory
	k = [get_headinfo(filelist[0]),get_headinfo(filelist[numfiles/2])]
	
	if verbose: print me+"Data extraction and calculations:",round(time.time()-t0,2),"seconds."
	
	##-------------------------------------------------------------------------
	## Plotting
	
	colour = ["b","r"]
	D_th = np.linspace(np.min(Delta),np.max(Delta),len(Delta)*20)
	
	fnt = 6
	fitxp = [0,0]
	fig, axs = plt.subplots(3,2, sharex=True)
	
	## Loop in Hop/not
	for i in [0,1]:
	
		ax = axs[0,0]
		ax.plot(Delta[i], S_fin[i], colour[i]+"o")
		ax.plot(D_th, SSS_theo(D_th**(2-i)), colour[i]+"--",\
			label="Optimal")
		fit = fit_SS(SSS_theo, Delta[i], S_fin[i]); fitxp[i]=round(fit[2],1)
		ax.plot(fit[0],fit[1], colour[i]+":", label="Fit ("+str(fitxp[i])+")")
		# ax.legend(prop={'size':fnt})
		ax.set_ylabel("$\Delta S_{\mathrm{SS}} / N\ln2$")
		ax.grid(i)
		
		ax = axs[0,1]
		ax.plot(Delta[i], t_SS[i], colour[i]+"o")
		ax.plot(Delta[i], N*t_SS_th[i], colour[i]+"--")
		ax.set_ylabel("$t_{\mathrm{SS}}$")
		ax.grid(i)
		ax.yaxis.major.formatter.set_powerlimits((0,0)) 
			
		ax = axs[1,0]
		ax.plot(Delta[i,1:], W_srt[i,1:], colour[i]+"o")
		ax.set_ylabel("$W_{\mathrm{total}}$ for sorting")
		ax.grid(i)
		ax.yaxis.major.formatter.set_powerlimits((0,0)) 
		
		ax = axs[2,0]
		ax.plot(Delta[i], Wdot_SS[i], colour[i]+"o")
		ax.plot(D_th, -SSW_theo(D_th,k[i],2-i), colour[i]+"--",\
			label="Optimal")
		# fit = fit_SS(SSW_theo, Delta[i], S_fin[i], k[i])
		# ax.plot(fit[0],fit[1], colour[i]+":",\
			# label="Fit: "+str(round(fit[2],1)))
		# ax.legend(loc="best",prop={'size':fnt})
		ax.set_xlim(left=1.0)
		ax.set_xlabel("$\Delta$")	
		ax.set_ylabel("$\dot W_{\mathrm{SS}}$")
		ax.grid(i)
		ax.yaxis.major.formatter.set_powerlimits((0,0))
	
	ax = axs[1,1]
	annotext = "Hopfield: BLUE\nNotfield: RED\nData: o\nTheory: --\nFit: .."
	ax.annotate(annotext,xy=(0.65,0.425),xycoords="figure fraction")
	plt.setp(ax.get_xticklabels(), visible=False)
	plt.setp(ax.get_yticklabels(), visible=False)
	
	
	ax = axs[2,1]
	errorplot(argv[1],ax,fitxp)
	ax.set_ylim(top=1.0)
	ax.set_xlim(left=1.0)
	ax.set_xlabel("$\Delta$")
		
	fig.suptitle("Hopfield and Notfield Properties Versus $\\Delta$")
	plt.tight_layout()
	plt.subplots_adjust(top=0.9)
	
	plt.savefig(plotfile)
	print me+"Plot saved to",plotfile
	
	print me+"Execution",round(time.time()-t0,2),"seconds."
	
	return


##=============================================================================
	
def fit_SS(func,x,y):
	"""
	Make a power-law fit to points y(x).
	Returns new x and y coordinates on finer grid.
	"""
	fitfunc = lambda x,nu: func(x,nu)
	popt, pcov = curve_fit(fitfunc, x, y)
	X = np.linspace(np.min(x),np.max(x),5*x.size)
	return X, fitfunc(X, *popt), popt[0]


##=============================================================================
##=============================================================================
if __name__=="__main__":
	main()
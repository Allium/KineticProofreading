
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
	ERR = np.zeros(numfiles)
	
	## Get data from all files
	for i in range(numfiles):
	
		Delta[i] = get_pars(filelist[i])
		k = get_headinfo(filelist[i])
	
		data = np.loadtxt(filelist[i], skiprows=10, unpack=True)
		data = data[:, ::int(data.shape[1]/npts)]
		
		t, ent, work, trans_C, trans_I, err = data[[0,5,6,8,9,12]]
		N = int(data[[1,2,3,4],0].sum())
		del data
		
		ent /= N*np.log(2)	
		SSidx = flatline(ent)
		
		## Assuming entropy is flat and work is linear
		
		S_fin[i] = np.mean(ent[SSidx:])
		t_SS[i] = t[SSidx]
		t_SS_th[i] = SSt_theo(k)
		# Wdot_srt[i] = np.mean(work[:SSidx-int(npts/20)])/t[SSidx]
		W_srt[i] = work[SSidx]
		Wdot_SS[i] = np.mean(work[SSidx:]-work[SSidx])/(t[-1]-t[SSidx])
		ERR[i] = err.mean()
	
	## ----------------------------------------------------
	
	## Separate into H and N
	newshape = (2,numfiles/2)
	Delta = Delta.reshape(newshape)
	S_fin = S_fin.reshape(newshape)
	t_SS = t_SS.reshape(newshape)
	t_SS_th = t_SS_th.reshape(newshape)
	W_srt = W_srt.reshape(newshape)
	Wdot_SS  = Wdot_SS.reshape(newshape)
	ERR = ERR.reshape(newshape)
	
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
	ERR = ERR[:,sortind[0]]

	## ----------------------------------------------------
	
	## Get k values, assuming the same for all files in directory
	k = [get_headinfo(filelist[0]),get_headinfo(filelist[numfiles/2])]
	
	if verbose: print me+"Data extraction and calculations:",round(time.time()-t0,2),"seconds."
	
	##-------------------------------------------------------------------------
	## Plotting
	
	fsl = 10
	colour = ["b","r","m"]
	label = ["Hopfield","Notfield"]
	D_th = np.linspace(np.min(Delta),np.max(Delta),len(Delta)*20)
	
	## SORTING ERROR RATE RATIO
	plt.figure("ERR"); ax = plt.gca()
	plotfile = argv[1]+"/DeltaPlot_0_ERR.png"
	plt.subplot(111)
	for i in [0,1]:
		fit = fit_par(ERR_fit, Delta[i], ERR[i])
		ax.plot(Delta[i], ERR[i], colour[i]+"o", label=label[i])
		ax.plot(Delta[i], Delta[i]**(fit[2]), colour[i]+":", label = "$\Delta^{%.2f}$" % fit[2])
		ax.plot(Delta[i], Delta[i]**(-2+i), colour[i]+"--", label = "$\Delta^{"+str(-2+i)+"}$")
	plt.xlim(left=1.0,right=Delta[0,-1])
	# plt.ylim(top=1.0)
	plt.xscale("log");	plt.yscale("log")
	plt.xlabel("$\\Delta$")
	plt.ylabel("Error Rate Ratio $\\langle\\dot I\\rangle/\\langle\\dot C\\rangle$")
	plt.grid()
	plt.legend(loc="upper right", fontsize=fsl)
	plt.savefig(plotfile)
	print me+"Plot saved to",plotfile
	
	## SS ENTROPY
	
	plt.figure("SSS"); ax = plt.gca()
	plotfile = argv[1]+"/DeltaPlot_1_SSS.png"
	for i in [0,1]:
		ax.plot(Delta[i], S_fin[i], colour[i]+"o", label=label[i])
		ax.plot(D_th, SSS_theo(D_th**(2-i)), colour[i]+"--",\
			label="Optimal")
		fit = fit_par(SSS_theo, Delta[i], S_fin[i])
		ax.plot(fit[0],fit[1], colour[i]+":", label="Fit ("+str(fit[2])+")")
	ax.set_xlim(left=1.0,right=Delta[0,-1])
	ax.set_xlabel("$\\Delta$")
	ax.set_ylabel("$\Delta S_{\mathrm{SS}} / N\ln2$")
	plt.grid()
	plt.legend(loc="upper right", fontsize=fsl)
	plt.savefig(plotfile); print me+"Plot saved to",plotfile
	
	## SS ENTROPY H/N
	plt.figure("SSSR"); ax = plt.gca()
	plotfile = argv[1]+"/DeltaPlot_2_SSSR.png"
	S_fin_ratio = (S_fin[1]+1.0)/(S_fin[0]+1.0)
	S_fin_th_ratio = (SSS_theo(Delta[1],1)+1.0)/(SSS_theo(Delta[0],2)+1.0)	
	ax.plot(Delta[0], S_fin_ratio, colour[2]+"o", label="Data")
	ax.plot(Delta[0], S_fin_th_ratio, colour[2]+"--",	label="Optimal")
	ax.set_xlim(left=1.0,right=Delta[0,-1])
	ax.set_ylim(bottom=0.0)
	ax.set_xlabel("$\\Delta$")
	ax.set_ylabel("$\\left(\\Delta S_{\\mathrm{SS}}^{\\mathrm{e}} + 1\\right)) /\
					\\left(\\Delta S_{\\mathrm{SS}}^{\\mathrm{p}} + 1\\right)$")
	plt.grid()
	plt.legend(loc="best")
	plt.savefig(plotfile)
	print me+"Plot saved to",plotfile
	
	## TIME TO REACH STEADY STATE
	
	plt.figure("tSS"); ax = plt.gca()
	plotfile = argv[1]+"/DeltaPlot_3_tSS.png"
	plt.subplot(111)
	for i in [0,1]:
		ax.plot(Delta[i], t_SS[i], colour[i]+"o", label=label[i])
		ax.plot(Delta[i,1:], N*t_SS_th[i,1:], colour[i]+"--", label="Theory")
	ax.set_xlim(left=1.0,right=Delta[0,-1])
	ax.set_xlabel("$\\Delta$")
	ax.set_ylabel("$t_{\mathrm{SS}}$")
	ax.yaxis.major.formatter.set_powerlimits((0,0)) 
	plt.grid()
	plt.legend(loc="best", fontsize=fsl)
	plt.savefig(plotfile); print me+"Plot saved to",plotfile
	
	## TOTAL WORK TO SORT
	
	plt.figure("Wsort"); ax = plt.gca()
	plotfile = argv[1]+"/DeltaPlot_4_Wsort.png"
	plt.subplot(111)
	for i in [0,1]:
		ax.plot(Delta[i,1:], W_srt[i,1:], colour[i]+"o", label=label[i])
	ax.set_xlim(left=1.0,right=Delta[0,-1])
	ax.set_xlabel("$\\Delta$")
	ax.set_ylabel("$W_{\mathrm{total}}$ for sorting")
	ax.yaxis.major.formatter.set_powerlimits((0,0)) 
	plt.grid()
	plt.legend(loc="lower right", fontsize=fsl)
	plt.savefig(plotfile); print me+"Plot saved to",plotfile
	
	## SS WORK RATE
	
	plt.figure("Wdot"); ax = plt.gca()
	plotfile = argv[1]+"/DeltaPlot_5_Wdot.png"
	plt.subplot(111)
	for i in [0,1]:
		ax.plot(Delta[i], Wdot_SS[i], colour[i]+"o", label=label[i])
		ax.plot(D_th, -SSW_theo(D_th,k[i],2-i), colour[i]+"--",\
			label="Theory")
	ax.set_xlim(left=1.0,right=Delta[0,-1])
	ax.set_xlabel("$\\Delta$")
	ax.set_ylabel("$\dot W_{\mathrm{SS}}$")
	ax.yaxis.major.formatter.set_powerlimits((0,0)) 
	plt.grid()
	plt.legend(loc="lower right", fontsize=fsl)
	plt.savefig(plotfile); print me+"Plot saved to",plotfile
		
	##	
	
	print me+"Execution",round(time.time()-t0,2),"seconds."
	
	return


##=============================================================================
	
def fit_par(func,x,y):
	"""
	Make a power-law fit to points y(x).
	Returns new x and y coordinates on finer grid.
	"""
	fitfunc = lambda x,nu: func(x,nu)
	popt, pcov = curve_fit(fitfunc, x, y)
	X = np.linspace(np.min(x),np.max(x),5*x.size)
	return X, fitfunc(X, *popt), popt[0]

def ERR_fit(D,nu):
	return D**nu

##=============================================================================
##=============================================================================
if __name__=="__main__":
	main()

import numpy as np
from matplotlib import pyplot as plt
from sys import argv
import os, glob, time
from SortTimePlot import get_pars, get_headinfo,\
	flatline, SSS_theo, SSW_theo

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
	
	Delta = np.zeros(numfiles)
	S_fin = np.zeros(numfiles)
	t_SS = np.zeros(numfiles)
	Wdot_srt = np.zeros(numfiles)
	Wdot_SS = np.zeros(numfiles)
	
	## Get data from all files
	for i in range(numfiles):
	
		Delta[i] = get_pars(filelist[i])
	
		data = np.loadtxt(filelist[i], skiprows=10, unpack=True)
		data = data[:, ::int(data.shape[1]/npts)]
		
		t, ent, work, trans_C, trans_I = data[[0,5,6,8,9]]
		N = int(data[[1,2,3,4],0].sum())
		del data
		
		ent /= N*np.log(2)	
		Dent, SSidx = flatline(ent, N, Delta)
		
		## Assuming entropy is flat and work is linear
		
		S_fin[i] = np.mean(ent[SSidx:])
		t_SS[i] = t[SSidx]
		Wdot_srt[i] = np.mean(work[:SSidx-int(npts/20)])/t[SSidx]
		Wdot_SS[i] = np.mean(work[SSidx:]-work[SSidx])/(t[-1]-t[SSidx])
	
	## Sort by Delta and reshape to separate Hop and Not
	sortind = np.argsort(Delta)
	newshape = [numfiles/2,2]
	Delta = Delta[sortind].reshape(newshape).T
	try:
		assert np.all(Delta[0]==Delta[1])
	except AssertionError:
		raise IOError(me+"Check files.\n"+filelist.tostring())
	S_fin = S_fin[sortind].reshape(newshape).T
	t_SS = t_SS[sortind].reshape(newshape).T
	Wdot_srt = Wdot_srt[sortind].reshape(newshape).T
	Wdot_SS = Wdot_SS[sortind].reshape(newshape).T
	
	## Get k values, assuming the same for all files in directory
	k = [get_headinfo(filelist[0]),get_headinfo(filelist[numfiles/2])]
	
	if verbose: print me+"Data extraction and calculations:",round(time.time()-t0,2),"seconds."
	
	##-------------------------------------------------------------------------
	## Plotting
	
	colour = ["b","r"]
	D_th = np.linspace(Delta[0,0],Delta[0,-1],len(Delta)*20)
	
	fig, axs = plt.subplots(3,2, sharex=True)
	
	for i in [0,1]:
		ax = axs[0,0]
		ax.plot(Delta[i], S_fin[i], colour[i]+"o")
		ax.plot(D_th, SSS_theo(D_th**(2-i)), colour[i]+"--",\
			label="Theory")
		ax.set_ylabel("$\Delta S_{\mathrm{SS}} / N\ln2$")
		ax.grid()
		
		ax = axs[0,1]
		ax.plot(Delta[i], t_SS[i], colour[i]+"o")	
		ax.set_ylabel("$t_{\mathrm{SS}}$")
		ax.grid()
		ax.yaxis.major.formatter.set_powerlimits((0,0)) 
			
		ax = axs[1,0]
		ax.plot(Delta[i], Wdot_srt[i], colour[i]+"o")
		ax.set_ylabel("$\dot W_{\mathrm{sort}}$")
		ax.grid()
		ax.yaxis.major.formatter.set_powerlimits((0,0)) 
		
		ax = axs[1,1]
		ax.plot(Delta[i], t_SS[i]*Wdot_srt[i], colour[i]+"o")
		ax.set_ylabel("$W_{\mathrm{total}}$ for sorting")
		ax.grid()
		ax.yaxis.major.formatter.set_powerlimits((0,0)) 
		
		ax = axs[2,0]
		ax.plot(Delta[i], Wdot_SS[i], colour[i]+"o")
		ax.plot(D_th, -SSW_theo(k[i],D_th,i), colour[i]+"--",\
			label="Theory")
		ax.set_xlabel("$\Delta$")	
		ax.set_ylabel("$\dot W_{\mathrm{SS}}$")
		ax.grid()
		ax.yaxis.major.formatter.set_powerlimits((0,0))
		
		ax = axs[2,1]
		ax.plot(Delta[i], np.zeros(len(Delta[i])), ".")
		ax.set_xlabel("$\Delta$")	
	
	# fig.suptitle("Hopfield and Notfield")
	plt.tight_layout()
	
	plt.savefig(plotfile)
	print me+"Plot saved to",plotfile
	
	print me+"Execution",round(time.time()-t0,2),"seconds."
	
	return



##=============================================================================
##=============================================================================
if __name__=="__main__":
	main()
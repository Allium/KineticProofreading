
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
	
	HISTORY
		2016/01/28	Started
	"""
	me = "SortDeltaPlot.main: "
	t0 = time.time()
	
	try: argv[1]
	except IndexError:
		print main.__doc__
		raise IndexError(me+"Please read docstring")
	
	plotfile = argv[1]+"/HopfieldDeltaPlots.png"

	##-------------------------------------------------------------------------

	filelist = glob.glob(argv[1]+"/*.txt")
	numfiles = len(filelist)
	
	Delta = np.zeros(numfiles)
	S_fin = np.zeros(numfiles)
	t_SS = np.zeros(numfiles)
	Wdot_srt = np.zeros(numfiles)
	Wdot_SS = np.zeros(numfiles)
	
	for i, datafile in enumerate(filelist):
	
		Delta[i] = get_pars(datafile)
	
		npts = 500
		data = np.loadtxt(datafile, skiprows=10, unpack=True)
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
	
	## Sort by Delta
	sortind = np.argsort(Delta)
	Delta = Delta[sortind]
	S_fin = S_fin[sortind]
	t_SS = t_SS[sortind]
	Wdot_srt = Wdot_srt[sortind]
	Wdot_SS = Wdot_SS[sortind]
	
	## Get k values, assuming the same for all files in directory
	k = get_headinfo(datafile)
	
	##-------------------------------------------------------------------------
	## Plotting
	
	fig, axs = plt.subplots(3,2, sharex=True)
	
	ax = axs[0,0]
	ax.plot(Delta, S_fin, "o")
	ax.plot(Delta, SSS_theo(Delta), color=ax.lines[-1].get_color(), ls="--",\
		label="Theory")
	# ax.set_xlabel("$\Delta$")
	ax.set_ylabel("$\Delta S_{\mathrm{SS}} / N\ln2$")
	ax.grid()
	
	ax = axs[0,1]
	ax.plot(Delta, t_SS, "o")	
	# ax.set_xlabel("$\Delta$")
	ax.set_ylabel("$t_{\mathrm{SS}}$")
	ax.grid()
	ax.yaxis.major.formatter.set_powerlimits((0,0)) 
		
	ax = axs[1,0]
	ax.plot(Delta, Wdot_srt, "o")
	# ax.set_xlabel("$\Delta$")	
	ax.set_ylabel("$\dot W_{\mathrm{sort}}$")
	ax.grid()
	ax.yaxis.major.formatter.set_powerlimits((0,0)) 
	
	ax = axs[1,1]
	ax.plot(Delta, t_SS*Wdot_srt, "o")
	# ax.set_xlabel("$\Delta$")
	ax.set_ylabel("$W_{\mathrm{total}}$ for sorting")
	ax.grid()
	ax.yaxis.major.formatter.set_powerlimits((0,0)) 
	
	ax = axs[2,0]
	ax.plot(Delta, Wdot_SS, "o")
	ax.plot(Delta, -SSW_theo(k,Delta), color=ax.lines[-1].get_color(), ls="--",\
		label="Theory")
	ax.set_xlabel("$\Delta$")	
	ax.set_ylabel("$\dot W_{\mathrm{SS}}$")
	ax.grid()
	ax.yaxis.major.formatter.set_powerlimits((0,0))
	
	ax = axs[2,1]
	ax.plot(Delta, np.zeros(len(Delta)), ".")
	ax.set_xlabel("$\Delta$")	
	
	fig.suptitle("Hopfield")
	plt.tight_layout()
	
	plt.savefig(plotfile)
	print me+"Plot saved to",plotfile
	
	print me+"Execution",round(time.time()-t0,2),"seconds."
	
	return



##=============================================================================
##=============================================================================
if __name__=="__main__":
	main()
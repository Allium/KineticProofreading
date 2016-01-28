
import numpy as np
from matplotlib import pyplot as plt
from sys import argv
import os, glob, time
from SortPlot import flatline, findpars

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
	except IndexError: print main.__doc__
	
	plotfile = argv[1]+"/HopfieldDeltaPlots.png"
	filelist = glob.glob(argv[1]+"/*.txt")
	numfiles = len(filelist)
	
	Delta = np.zeros(numfiles)
	S_fin = np.zeros(numfiles)
	Wdot_evo = np.zeros(numfiles)
	Wdot_SS = np.zeros(numfiles)
	t_SS = np.zeros(numfiles)
	
	for i, datafile in enumerate(filelist):
	
		Delta[i] = findpars(datafile)
	
		npts = 500
		data = np.loadtxt(datafile, skiprows=4, unpack=True)
		data = data[:, ::int(data.shape[1]/npts)]
		
		t, ent, work, tran = data
		Dent, SSidx = flatline(ent, Delta[i])
		
		## Assuming entropy is flat and work is linear
		
		S_fin[i] = np.mean(ent[SSidx:])
		Wdot_evo[i] = np.mean(work[:SSidx-int(npts/20)])/t[SSidx]
		Wdot_SS[i] = np.mean(work[SSidx:]-work[SSidx])/(t[-1]-t[SSidx])
		t_SS[i] = t[SSidx]

	sortind = np.argsort(Delta)
	Delta = Delta[sortind]
	S_fin = S_fin[sortind]
	Wdot_evo = Wdot_evo[sortind]
	Wdot_SS = Wdot_SS[sortind]
	t_SS = t_SS[sortind]
	
	fig, ax = plt.subplots(2,2)
	ax[0,0].plot(Delta, S_fin, "o")
	ax[0,0].set_xlabel("$\Delta$")
	ax[0,0].set_ylabel("$S_{final}$")
	ax[0,0].grid()
	
	ax[0,1].plot(Delta, Wdot_evo, "o")
	ax[0,1].set_xlabel("$\Delta$")	
	ax[0,1].set_ylabel("$\dot W_{evo}$")
	ax[0,1].grid()
	ax[0,1].yaxis.major.formatter.set_powerlimits((0,0)) 
	
	ax[1,0].plot(Delta, Wdot_SS, "o")
	ax[1,0].set_xlabel("$\Delta$")	
	ax[1,0].set_ylabel("$\dot W_{SS}$")
	ax[1,0].grid()
	ax[1,0].yaxis.major.formatter.set_powerlimits((0,0))
	
	ax[1,1].plot(Delta, t_SS,"o")	
	ax[1,1].set_xlabel("$\Delta$")
	ax[1,1].set_ylabel("$t_{SS}$")
	ax[1,1].grid()
	ax[1,1].yaxis.major.formatter.set_powerlimits((0,0)) 
	
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

import numpy as np
from scipy.signal import fftconvolve
from scipy.optimize import curve_fit
from scipy.interpolate import Rbf
from matplotlib import pyplot as plt
import os, glob, optparse, time
from SymPlot import *

##=============================================================================
def main():
	"""
	"""
	me = "PhaseDiagram.main: "
	t0 = time.time()

	## Options
	parser = optparse.OptionParser()
	parser.add_option('-D','--Delta',
		dest="Delta", default=5.0, type="float")
	parser.add_option('--plot_theory',
		dest="theory", default=0, type="int")
	parser.add_option('-s','--show',
		dest="showfig", default=False, action="store_true")
	parser.add_option('--nosave',
		dest="nosave", default=False, action="store_true")
	parser.add_option('-a','--plotall',
		dest="plotall", default=False, action="store_true")
	parser.add_option('-v','--verbose',
		dest="verbose", default=False, action="store_true")
	opt, args = parser.parse_args()
	
	Delta = opt.Delta
	theory = opt.theory
	showfig = opt.showfig
	nosave = opt.nosave
	plotall = opt.plotall
	verbose = opt.verbose

	assert os.path.isdir(args[0]), me+"First CLA must be a directory path."
	assert theory == 0 or 1 or 2, me+"plot_theory = 0,1,2."
	
	if plotall:
		showfig = False
		raise AttributeError, me+"Functionality doesn't exist."
	else:
		phase_delta(args[0], Delta, theory, nosave, verbose)
	
	if verbose: print me+"Execution",round(time.time()-t0,2),"seconds."
	if showfig: plt.show()
	
	return
	

##=============================================================================
def phase_delta(dirpath, Delta, theory, ns, vb):
	
	me = "PhaseDiagram.phase_delta: "
	t0 = time.time()

	##-------------------------------------------------------------------------
	
	## Select hopfield files
	filelist = np.sort(glob.glob(dirpath+"/*Hop*_"+str(Delta)+"_*.txt"))
	numfiles = len(filelist)
	if vb: print me+"Found",numfiles,"files with Delta =",Delta

	## Initialise arrays
	DS	= np.zeros([2,numfiles])
	T	= np.zeros([2,numfiles])
	W_srt	= np.zeros([2,numfiles])
	Wdot	= np.zeros([2,numfiles])
	
	## Get data from all files
	for i in range(numfiles):
	
		k = get_headinfo(filelist[i])
		data = np.loadtxt(filelist[i], skiprows=10, unpack=True)
		
		## Read relevant columns
		N = int(data[[1,2,3,4],0].sum())
		t, A1, work = data[[0,1,6]]
		
		## Normalise ent and work
		ent = calc_ent_norm(A1/(N/2.0))
		work *= (0.0 if Delta == 1.0 else 1.0/(np.log(Delta) * N*np.log(2)))
		
		## SS index
		SSidx = (0 if Delta == 1.0 else flatline(ent))
		tSS = t[SSidx]
		
		## Simulation calculations
		if theory == 0 or theory == 1:
			DS[0,i]	= np.mean(ent[SSidx:])
			T[0,i]	= tSS
			Wdot[0,i]	= np.mean(work[SSidx:]-work[SSidx])/(t[-1]-tSS)
			W_srt[0,i]	= work[SSidx]
		
		## Theory calculations
		if theory == 1 or theory == 2:
			DS[1,i]	= SSS_theo(Delta,k)
			T[1,i]	= SSt_theo(Delta,k)*N
			Wdot[1,i]	= -SSW_theo(Delta,k)/(N*np.log(2))
		
	## Sort
	idx = np.argsort(T[0])
	T = T[:,idx]
	DS = DS[:,idx]
	Wdot = Wdot[:,idx]
	W_srt = W_srt[:,idx]
	
	if vb: print me+"Data extraction and calculations:",round(time.time()-t0,2),"seconds."
	
	##-------------------------------------------------------------------------
	## Plotting
	
	fsa,fsl,fst = 16,14,18
	
	plotkey = ["sim","simpre","pre"]
	plotfile = dirpath+"/PhaseDiagram_"+plotkey[theory]+"_D"+str(Delta)+".png"
	
	plt.figure("Phase Diagram")
	plt.subplot(111)
	
	ent = True
	if ent:
		z = DS
		tit = "$\\Delta S$. "
	else:
		z = DS - W_srt
		tit = "$\\Delta S - W_{\\rm sort}$"
	
	"""## Grid of interpolation points
	T_g, Wdot_g = np.meshgrid(np.linspace(T.min(), T.max(), 10), np.linspace(Wdot.min(), Wdot.max(), 10))
	## Interpolate
	z_interpfunc = Rbf(T, Wdot, z, function="linear")
	z_g = z_interpfunc(T_g, Wdot_g)
	plt.imshow(z_g, vmin=DS.min(), vmax=DS.max(), origin="lower", extent=[Wdot.min(),Wdot.max(),T.min(),T.max()])"""

	if theory == 0 or theory == 1:
		plt.scatter(T[0], Wdot[0], c=z[0], marker="o", label="Simulation")
	if theory == 1 or theory == 2:
		plt.scatter(T[1], Wdot[1], c=z[1], marker="x", label="Theory")
	
	plt.xlim(left=0.0)
	plt.ylim(top=0.0,bottom=Wdot[0].min())
	plt.colorbar()
	if ent:	plt.clim(0.0,-1.0)
	
	plt.xlabel("$\\tau$",fontsize=fsa)
	plt.ylabel("$\\dot W_{\\rm SS}$",fontsize=fsa)
	plt.suptitle(tit+"$\\Delta = "+str(Delta)+"$.",fontsize=fst)
	plt.legend(loc="lower right",fontsize=fsl)
	plt.gca().xaxis.major.formatter.set_powerlimits((0,0)) 
	plt.gca().yaxis.major.formatter.set_powerlimits((0,0)) 
	plt.grid()

	if not ns:
		plt.savefig(plotfile)
		if vb: print me+"Plot saved to",plotfile

	if vb: print me+"Execution",round(time.time()-t0,2),"seconds."
	
	return
	


##=============================================================================
##=============================================================================
if __name__=="__main__":
	main()
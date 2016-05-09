
import numpy as np
from scipy.signal import fftconvolve
from scipy.optimize import curve_fit
from scipy.interpolate import Rbf
from matplotlib import pyplot as plt
import os, glob, optparse, time, sys
from SymPlot import *

##=============================================================================
def main():
	"""
	NAME
		PhaseDiagram.py
		
	EXECUTION
		python PhaseDiagram.py dirpath opts/flags
		
	OPTS
			--plot_theory	0	0,1,2
	
	FLAGS
		-a	--plotall	False
		-s	--show		False
			--nosave	False
		-v	--verbose	False
		-h	--help		False
	"""
	me = "PhaseDiagram.main: "
	t0 = time.time()

	## Options
	parser = optparse.OptionParser(conflict_handler="resolve")
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
	parser.add_option('-h','--help',
		dest="help", default=False, action="store_true")
	opt, args = parser.parse_args()
	
	if opt.help:
		print main.__doc__
		return	
	theory = opt.theory
	showfig = opt.showfig
	nosave = opt.nosave
	plotall = opt.plotall
	verbose = opt.verbose

	args[0] = args[0].replace("\\","/")
	assert os.path.isdir(args[0]), me+"First CLA must be a directory path."
	assert theory == 0 or 1 or 2, me+"plot_theory = 0,1,2."
	
	try:				Delta = float(args[0][args[0].find("/D")+2:])
	except ValueError:	Delta = float(args[0][args[0].find("/D")+2:-1])
	
	if plotall:
		showfig = False
		raise AttributeError, me+"Functionality doesn't exist."
	elif theory != 2:
		phase_dat = phase_calc(args[0], Delta, theory, verbose)
		phase_plot(args[0], theory, nosave, verbose, Delta, *phase_dat)
	elif theory == 2:
		phase_dat = phase_theo(Delta, verbose)
		phase_plot(args[0], theory, nosave, verbose, Delta, *phase_dat)
	
	if verbose: print me+"Execution",round(time.time()-t0,2),"seconds."
	if showfig: plt.show()
	
	return
	

##=============================================================================
def phase_calc(dirpath, Delta, th, vb):
	
	me = "PhaseDiagram.phase_calc: "
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
	j=0
	
	## Get data from all files
	for i in range(numfiles):
	
		if vb and i%(numfiles/20)==0:
			sys.stdout.write("\r"+me+"Processing: [%-20s] %d%%" % ('='*j, 5*j))
			j+=1
	
		k = get_headinfo(filelist[i])
		data = np.loadtxt(filelist[i], skiprows=10, unpack=True)
		
		## Read relevant columns
		N = int(data[[1,2,3,4],0].sum())
		t, A1, work = data[[0,1,6]]
		
		## Normalise time, ent and work
		t /= N
		ent = calc_ent_norm(A1/(N/2.0))
		work *= (0.0 if Delta == 1.0 else 1.0/(np.log(Delta) * N*np.log(2)))
		
		## SS index
		SSidx = (0 if Delta == 1.0 else flatline(ent))
		tSS = t[SSidx]
		
		## Simulation calculations
		DS[0,i]	= np.mean(ent[SSidx:])
		T[0,i]	= tSS
		Wdot[0,i]	= np.mean(work[SSidx:]-work[SSidx])/(t[-1]-tSS)
		W_srt[0,i]	= work[SSidx]
		
		## Theory calculations
		if th == 1:
			DS[1,i]	= SSS_theo(Delta,k)
			T[1,i]	= SSt_theo(Delta,k)
			Wdot[1,i]	= -SSW_theo(Delta,k)
		
	## Sort
	idx = np.argsort(T[0])
	T = T[:,idx]
	DS = DS[:,idx]
	Wdot = Wdot[:,idx]
	W_srt = W_srt[:,idx]
	
	if vb: print "\n"+me+"Data extraction and calculations:",round(time.time()-t0,2),"seconds."
	
	return (T, DS, Wdot, W_srt)

	
##=============================================================================	
def phase_plot(dirpath, th, ns, vb, Delta, T, DS, Wdot, W_srt):
	
	me = "PhaseDiagram.phase_plot: "
	t0 = time.time()
	
	fsa,fsl,fst = 16,14,18
	
	## Z-AXIS
	
	ent = True
	if ent:
		z = DS
		tit = "$\\Delta S$. "
		plotkey = "_S"
	else:
		## Remember W<0
		if th != 0: print me+"Warning: No theory for W_sort."
		z = DS - W_srt
		tit = "$\\Delta S - W_{\\rm sort}$. "
		plotkey = "_SW"
	plotkey = ["sim","simpre","pre"][th] + plotkey
	plotfile = dirpath+"/PhaseDiagram_"+plotkey+"_D"+str(Delta)+".png"
	
	## PLOTTING
	
	plt.figure("Phase Diagram")
	plt.subplot(111)
	
	"""## Grid of interpolation points
	T_g, Wdot_g = np.meshgrid(np.linspace(T.min(), T.max(), 10), np.linspace(Wdot.min(), Wdot.max(), 10))
	## Interpolate
	z_interpfunc = Rbf(T, Wdot, z, function="linear")
	z_g = z_interpfunc(T_g, Wdot_g)
	plt.imshow(z_g, vmin=DS.min(), vmax=DS.max(), origin="lower", extent=[Wdot.min(),Wdot.max(),T.min(),T.max()])"""

	if th == 0 or th == 1:
		plt.scatter(T[0], Wdot[0], c=z[0], marker="o", label="Simulation")
	if th == 1 or th == 2:
		plt.scatter(T[1], Wdot[1], c=z[1], marker="x", label="Theory")
	
	## PLOT PROPERTIES
	
	plt.xlim(left=0.0)
	plt.ylim(top=0.0,bottom=Wdot[0].min())
	plt.colorbar()
	plt.clim(vmin=-1.0)
	if ent:	plt.clim(vmax=0.0)
	
	plt.xlabel("$\\tau$",fontsize=fsa)
	plt.ylabel("$\\dot W_{\\rm SS}$",fontsize=fsa)
	plt.suptitle(tit+"$\\Delta = "+str(Delta)+"$.",fontsize=fst)
	plt.legend(loc="lower right",fontsize=fsl)
	plt.gca().xaxis.major.formatter.set_powerlimits((0,0)) 
	plt.gca().yaxis.major.formatter.set_powerlimits((0,0)) 
	plt.grid()

	## SAVING
	
	if not ns:
		plt.savefig(plotfile)
		if vb: print me+"Plot saved to",plotfile

	if vb: print me+"Plotting",round(time.time()-t0,2),"seconds."
	
	return
	

##=============================================================================
def phase_theo(Delta, vb):

	me = "PhaseDiagram.phase_theo: "
	t0 = time.time()
	
	print me+"Plotting theoretical phase diagram for Delta =",Delta
	
	A1B = np.arange(0.005,0.04,0.005)
	BA1 = np.arange(0.001,0.02,0.002)
	BC  = np.arange(0.001,0.04,0.002)
	CA1 = np.arange(0.001,0.02,0.002)
	klist = [{"A1B1": a1b, "B1A1": ba1, "B1C1": bc, "C1A1": ca1} \
				for a1b in A1B for ba1 in BA1 for bc in BC for ca1 in CA1]
	npoints = len(klist)
	if vb:	print me+"Computing",npoints,"points."

	## Initialise arrays
	DS	= np.zeros([2,npoints])
	T	= np.zeros([2,npoints])
	W_srt	= np.zeros([2,npoints])
	Wdot	= np.zeros([2,npoints])
	
	for i,k in enumerate(klist):
		DS[1,i]	= SSS_theo(Delta,k)
		T[1,i]	= SSt_theo(Delta,k)
		Wdot[1,i]	= -SSW_theo(Delta,k)
		
	## Sort
	idx = np.argsort(T[1])
	T = T[:,idx]
	DS = DS[:,idx]
	Wdot = Wdot[:,idx]
	
	Wdot[0,0] = Wdot.min()
	
	if vb: print me+"Calculations:",round(time.time()-t0,2),"seconds."
	
	return (T, DS, Wdot, W_srt)

##=============================================================================
##=============================================================================
if __name__=="__main__":
	main()
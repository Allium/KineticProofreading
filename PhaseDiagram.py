
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
		-S	--entropy	False	Use entropy as quality metric (default S-W).
		-s	--show		False
			--nosave	False
		-a	--allplot	False	Make all relevant plots for directory.
		-v	--verbose	False
		-h	--help		False
	"""
	
	me = "PhaseDiagram.main: "
	t0 = time.time()

	## Options
	parser = optparse.OptionParser(conflict_handler="resolve")
	parser.add_option('--plot_theory',
		dest="theory", default=0, type="int")
	parser.add_option('-S','--entropy',
		dest="plotent", default=False, action="store_true")
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
	plotent = opt.plotent
	theory 	= opt.theory
	showfig = opt.showfig
	nosave 	= opt.nosave
	plotall	= opt.plotall
	verbose = opt.verbose
	dirpath = args[0]

	## Tests
	dirpath = dirpath.replace("\\","/")
	assert os.path.isdir(dirpath), me+"First CLA must be a directory path."
	assert theory == 0 or 1 or 2, me+"Usage: --plot_theory [0,1,2]."
	if (plotent is False and theory != 0):
		print me+"Warning! No theory for W_sort. Plotting entropy change only."
		plotent = True
	
	## Infer Delta from dirpath
	try:				Delta = float(dirpath[dirpath.find("/D")+2:])
	except ValueError:	Delta = float(dirpath[dirpath.find("/D")+2:-1])
	
	## Calculations
	if plotall:
		if verbose: print me+"Plotting all phase diagrams for Delta =",Delta,"\n"
		plot_all(dirpath, Delta, verbose)
		return
	elif theory == 0 or theory == 1:
		if verbose: print me+"Plotting simulated phase diagram for Delta =",Delta
		phase_dat = phase_calc(dirpath, Delta, theory, verbose)
	elif theory == 2:
		if verbose: print me+"Plotting theoretical phase diagram for Delta =",Delta
		phase_dat = phase_theo(Delta, verbose)
	
	## Plotting
	phase_plot(dirpath, theory, plotent, nosave, verbose, Delta, *phase_dat)
	
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
	
		if vb and i%((numfiles-1)/20)==0:
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
		
		if (tSS>3e2 and Wdot[0,i]<-0.8e-2): print [i],filelist[i],k
		
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
def phase_plot(dirpath, th, pe, ns, vb, Delta, T, DS, Wdot, W_srt):
	
	me = "PhaseDiagram.phase_plot: "
	t0 = time.time()
	
	fsa,fsl,fst = 16,14,18
	
	## Z-AXIS
	
	if pe:
		z = DS
		tit = "$\\Delta S$. "
		plotkey = "_S"
	else:
		## Remember W<0
		z = DS - W_srt
		tit = "$\\Delta S - W_{\\rm sort}$. "
		plotkey = "_SW"
	plotkey = ["sim","simpre","pre"][th] + plotkey
	plotfile = dirpath+"/PhaseDiagram_"+plotkey+"_D"+str(Delta)+".png"
	
	## PLOTTING
	
	fig = plt.figure("Phase Diagram")
	plt.subplot(111)
	
	if th != 0:	plt.scatter(T[1], Wdot[1], marker="x", c=z[1], label="Theory")
	plt.scatter(T[0], Wdot[0], marker="o", c=z[0], edgecolor="None", label="Simulation")
	
	## PLOT PROPERTIES
	
	plt.xlim(left=0.0)
	plt.ylim(top=0.0,bottom=Wdot[0].min())
	plt.colorbar()
	plt.clim(vmin=-1.0)
	if pe:	plt.clim(vmax=0.0)
	
	plt.xlabel("$\\tau$",fontsize=fsa)
	plt.ylabel("$\\dot W_{\\rm SS}$",fontsize=fsa)
	plt.suptitle(tit+"$\\Delta = "+str(Delta)+"$.",fontsize=fst)
	plt.legend(loc="lower right",fontsize=fsl)
	plt.gca().xaxis.major.formatter.set_powerlimits((0,0)) 
	plt.gca().yaxis.major.formatter.set_powerlimits((0,0)) 
	plt.grid()

	## SAVING
	
	if not ns:
		plt.savefig(plotfile, dpi=2*fig.dpi)
		if vb: print me+"Plot saved to",plotfile

	if vb: print me+"Plotting",round(time.time()-t0,2),"seconds."
	
	return
	

##=============================================================================
def phase_theo(Delta, vb):

	me = "PhaseDiagram.phase_theo: "
	t0 = time.time()
	
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
def plot_all(dirpath,Delta,verbose):
	## S theory
	phase_dat = phase_theo(Delta, verbose)
	phase_plot(dirpath, 2, True, False, verbose, Delta, *phase_dat)
	plt.close();	print "\n"
	## S sim
	phase_dat = phase_calc(dirpath, Delta, 0, verbose)
	phase_plot(dirpath, 0, True, False, verbose, Delta, *phase_dat)
	plt.close();	print "\n"
	## S sim+theo
	phase_dat = phase_calc(dirpath, Delta, 1, verbose)
	phase_plot(dirpath, 1, True, False, verbose, Delta, *phase_dat)
	plt.close();	print "\n"
	## WS sim
	phase_dat = phase_calc(dirpath, Delta, 0, verbose)
	phase_plot(dirpath, 0, False, False, verbose, Delta, *phase_dat)
	return exit()

##=============================================================================
##=============================================================================
if __name__=="__main__":
	main()
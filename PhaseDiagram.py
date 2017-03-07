me0 = "PhaseDiagram"

import numpy as np
import os, glob, optparse, time, sys

import matplotlib as mpl
if "SSH_TTY" in os.environ:
	print me0+": Using Agg backend."
	mpl.use("Agg")
from matplotlib import pyplot as plt

from SymPlot import *

from Utils import fs, set_mplrc, filename_par
set_mplrc(fs)

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
	
	me = me0+".main: "
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
	# if (plotent is False and theory != 0):
		# print me+"Warning! No theory for W_sort. Plotting entropy change only."
		# plotent = True
	
	## Infer Delta from dirpath
	try:				Delta = filename_par(dirpath, "/D")#float(dirpath[dirpath.find("/D")+2:])
	except ValueError:	Delta = float(dirpath[dirpath.find("/D")+2:-1])
	
	## Calculations
	if plotall:
		if verbose: print me+"Plotting all phase diagrams for Delta =",Delta,"\n"
		plot_all(dirpath, Delta, showfig, verbose)
		return
	
	## Plotting
	phase_plot(dirpath, theory, plotent, nosave, verbose, Delta)
	
	if verbose: print me+"Execution",round(time.time()-t0,2),"seconds."
	if showfig: plt.show()
	
	return
	

##=============================================================================
def phase_calc(dirpath, Delta, th, vb):
	
	me = "PhaseDiagram.phase_calc: "
	t0 = time.time()

	##-------------------------------------------------------------------------
	
	if th==2:
		print me+"Transferring to theoretical calculation."
		return phase_theo(Delta, vb)
		
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
		try:	data = np.loadtxt(filelist[i], skiprows=10, unpack=True)
		except ValueError:	raise IOError, me+"Issue with file "+filelist[i]
		
		## Read relevant columns
		N = int(data[[1,2,3,4],0].sum())
		t, A1, work = data[[0,1,6]]
		
		## Normalise time, ent and work
		t /= N
		# t *= k["B1C1"]
		ent = calc_ent_norm(A1/(N/2.0))
		work *= (0.0 if Delta == 1.0 else (np.log(Delta) / N))
		
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
		
	## SAVING
	calcfile = dirpath+"/CALC_D"+str(int(Delta))+".npz"
	np.savez(calcfile, Delta=Delta, DS=DS, T=T, Wdot=Wdot, W_srt=W_srt)
	if vb:
		print me+"Calculations saved to",calcfile

	return {"Delta":Delta,"DS":DS,"T":T,"Wdot":Wdot,"W_srt":W_srt}
	
	
##=============================================================================
def phase_theo(Delta, vb):
	"""
	Calculate theoretical data.
	"""

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
		W_srt[1,i]	= Wsort_theo(Delta,k)
		Wdot[1,i]	= -SSW_theo(Delta,k)
		
	## Sort
	idx = np.argsort(T[1])
	T = T[:,idx]
	DS = DS[:,idx]
	Wdot = Wdot[:,idx]
	W_srt = W_srt[:,idx]
	
	Wdot[0,0] = Wdot.min()
	
	if vb: print me+"Calculations:",round(time.time()-t0,2),"seconds."
	
	return {"Delta":Delta,"T":T,"DS":DS,"Wdot":Wdot,"W_srt":W_srt}
	
	
##=============================================================================	
def phase_plot(dirpath, th, pe, nosave, vb, Delta):
	"""
	Arguments are:
	dirpath		path to directory
	th			0,2,1 plot only sim, theory, both
	pe			plot entropy rather than S-W
	"""
	me = "PhaseDiagram.phase_plot: "
	t0 = time.time()
	
	##-------------------------------------------------------------------------
	## Read in existing data or calculate afresh
		
	try:
		assert th==0
		data = np.load(dirpath+"/CALC_D"+str(int(Delta))+".npz")
		print me+"Data file found:",dirpath+"/CALC_D"+str(int(Delta))+".npz"
	except (IOError, AssertionError):
		print me+"No data found. Calculating."
		data = phase_calc(dirpath, Delta, th, vb)
		
	DS = data["DS"]
	T = data["T"]
	W_srt = data["W_srt"]
	Wdot = data["Wdot"]
	del data
	
	##-------------------------------------------------------------------------
	
	fs["saveext"]="jpg"
	
	## Z-AXIS
	
	if pe:
		z = DS
		tit = r"$\Delta S/N\ln2$. "
		plotkey = "_S"
	else:
		## Remember W<0
		z = DS - W_srt/np.log(2)
		tit = r"$(\Delta S - W^{\rm sort}/T)/N\ln2$. "
		plotkey = "_SW"
	plotkey = ["sim","simpre","pre"][th] + plotkey
	plotfile = dirpath+"/PhaseDiagram_"+plotkey+"_D"+str(int(Delta))+"."+fs["saveext"]
	
	## PLOTTING
	
	fig, ax = plt.subplots(1,1, figsize=fs["figsize"])
	
	## Plot Wsort vs DS with colour according to tau
	if 1:
		plotfile = dirpath+"/PhaseDiagram2_"+plotkey+"_D"+str(int(Delta))+"."+fs["saveext"]
		im = plt.scatter(DS[th/2], W_srt[th/2], marker="o", c=T[th/2], edgecolor="None", label="Simulation")
		Seq = SSS_theo(Delta, {"A1B1":0.1,"B1C1":0.1,"C1A1":0.0,"B1A1":0.1})
		plt.scatter(Seq,0.0, marker="o", c="k", s=40, label="Passive")
		ax.set_xlim(-1,0)
		ax.set_ylim(top=0.0)
		ax.set_ylim(bottom=-14.0)
		ax.set_xlabel(r"$\Delta S/N\ln2$")
		ax.set_ylabel(r"$W^{\rm sort}/NT$")
		cbar = fig.colorbar(im, ax=ax,aspect=50)
		cbar.ax.tick_params(labelsize=fs["fsl"]-4)
		cbar.locator = MaxNLocator(3);	cbar.update_ticks()
		plt.subplots_adjust(right=1.0)
	
	## Plot Wdot vs tau with colour according to DS OR DS-W
	else:
		if th != 0:	plt.scatter(T[1], Wdot[1], marker="x", c=z[1], label="Theory")
		if th != 2: plt.scatter(T[0], Wdot[0], marker="o", c=z[0], edgecolor="None", label="Simulation")
		
		# PLOT PROPERTIES
		
		ax.set_xlim(left=0.0)
		ax.set_ylim(top=0.0,bottom=Wdot[0].min())
		plt.colorbar(ax=ax)
		plt.clim(vmin=-1.0)
		if pe:	plt.clim(vmax=0.0)
		
		if int(Delta)==10:
			ax.set_xlim(0.0,6e2)
			ax.set_ylim(-2e-2,0.0)
		elif int(Delta)==5:
			ax.set_xlim(0.0,9e2)
			ax.set_ylim(-1.7e-2,0.0)
			
		ax.set_xlabel(r"$\tau/N$")
		ax.set_ylabel(r"$\dot W^{\rm SS}/NT$")
	
	
	ax.xaxis.major.formatter.set_powerlimits((0,0)) 
	ax.yaxis.major.formatter.set_powerlimits((0,0)) 
	ax.grid()
	plt.subplots_adjust(left=0.15)

	## SAVING
	
	if not nosave:
		fig.savefig(plotfile, dpi=2*fig.dpi)
		if vb: print me+"Plot saved to",plotfile

	if vb: print me+"Plotting",round(time.time()-t0,2),"seconds."
	
	return
	


##=============================================================================
def plot_all(dirpath,Delta,showfig,verbose):

#	## S theory
	phase_plot(dirpath, 2, True, False, verbose, Delta)
#	## WS theory
#	phase_plot(dirpath, 2, False, False, verbose, Delta)

	## S sim
	phase_plot(dirpath, 0, True, False, verbose, Delta)
	# WS sim
	# phase_plot(dirpath, 0, False, False, verbose, Delta)
	
	##
	if showfig: plt.show()
	return

##=============================================================================
##=============================================================================
if __name__=="__main__":
	main()

me0 = "Trade"

import numpy as np
from matplotlib import pyplot as plt
from cycler import cycler
import os, glob, optparse, time, sys

from PhaseDiagram import phase_calc

from Utils import fs, set_mplrc, filename_par
set_mplrc(fs)

##=============================================================================
def main():
	"""
	For a given Delta and quality (Delta S), what is trade-off between work and time?
	
	To do:
		Pass entropy values as CLAs
	"""
	
	me = me0+".main: "
	t0 = time.time()

	## Options
	parser = optparse.OptionParser(conflict_handler="resolve")
	parser.add_option('-s','--show',
		dest="showfig", default=False, action="store_true")
	parser.add_option('--nosave',
		dest="nosave", default=False, action="store_true")
	parser.add_option('-v','--verbose',
		dest="verbose", default=False, action="store_true")
	parser.add_option('-h','--help',
		dest="help", default=False, action="store_true")
	opt, args = parser.parse_args()
	
	if opt.help:
		print main.__doc__
		return
	showfig = opt.showfig
	nosave 	= opt.nosave
	verbose = opt.verbose
	dirpath = args[0]

	## Tests
	assert os.path.isdir(dirpath), me+"First CLA must be a directory path."
	dirpath = dirpath.replace("\\","/")
		
	## Plotting
	# entropy = [0.40, 0.56, 0.70]
	entropy = [0.2,0.3,0.4,0.5,0.6,0.7,0.8]
	trade_plot(dirpath, entropy, nosave, verbose)
	
	if verbose: print me+"Execution",round(time.time()-t0,2),"seconds."
	if showfig: plt.show()
	
	return
	
	
##=============================================================================	
def trade_plot(dirpath, entropy, nosave, vb):
	"""
	"""
	me = me0+".trade_plot: "
	t0 = time.time()
	
	## Infer Delta from dirpath
	Delta = filename_par(dirpath, "/D")
	
	##-------------------------------------------------------------------------
	## Read in existing data or calculate afresh using phase_calc
		
	try:
		data = np.load(dirpath+"/CALC_D"+str(int(Delta))+".npz")
		print me+"Data file found:",dirpath+"/CALC_D"+str(int(Delta))+".npz"
	except IOError:
		print me+"No data found. Calculating."
		data = phase_calc(dirpath, Delta, 0, vb)
		
	DS = data["DS"][0]
	T = data["T"][0]
	W_srt = data["W_srt"][0]
	Wdot = data["Wdot"][0]
	del data
	
	##-------------------------------------------------------------------------
	## ORDER DATA
		
	## Order data by value of DS
	idx = DS.argsort()
	DS	= DS[idx]
	T	= T[idx]
	W_srt = W_srt[idx]
	Wdot  = Wdot[idx]
	
	## Ensure entropy is negative
	entropy = [-abs(Si) for Si in entropy]
	
	##-------------------------------------------------------------------------
	## PARAMETERS
	
	## Entropy range
	tol = 0.02
	
	## Outplot
	fs["saveext"]="jpg"
	plotfile = dirpath+"/Trade_D"+str(int(Delta))+"."+fs["saveext"]
		
	##-------------------------------------------------------------------------
	## PLOT DATA

	fig, ax = plt.subplots(1,1, figsize=fs["figsize"])
	
	## Colours
	colours = plt.cm.jet(np.linspace(0,1,len(entropy)))[::-1]
	ax.set_prop_cycle(cycler('color', colours))
	
	for Si in entropy:
		idx = (DS>Si-tol)*(DS<Si+tol)
		ax.plot(T[idx], W_srt[idx]/np.log(2), "o", ms=6, label=r"$%.2f$"%(Si))
	
	##-------------------------------------------------------------------------
	## FIGURE
	
	ax.set_xlim(left=0.0)
	ax.set_ylim(top=0.0)
	ax.xaxis.major.formatter.set_powerlimits((0,0)) 
	ax.yaxis.major.formatter.set_powerlimits((0,0))
	
	ax.set_xlabel(r"$\tau\times bc/N$")
	# ax.set_xlabel(r"$\tau/N$")
	ax.set_ylabel(r"$W^{\rm sort}/NT\ln2$")
	
	leg = ax.legend(loc="best")#, ncol=min(3,len(entropy)/2))
	leg.set_title(r"$\Delta S/N\ln2$", prop={"size":fs["fsl"]})
	leg.get_frame().set_alpha(0.7)
	ax.grid()

	##-------------------------------------------------------------------------
	## SAVING
	
	if not nosave:
		fig.savefig(plotfile)
		if vb: print me+"Plot saved to",plotfile

	if vb: print me+"Plotting",round(time.time()-t0,2),"seconds."
	
	return
	
##=============================================================================
##=============================================================================
if __name__=="__main__":
	main()

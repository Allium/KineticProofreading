
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d
from matplotlib import pyplot as plt
from sys import argv
import os, glob, time

##=============================================================================
def main():
	"""
	NAME
		SortTimePlot.py
	
	PURPOSE
	
	EXECUTION
		python SortTimwPlot.py [filename]
		
	OUTPUT
		png plot saved in file directory
	
	EXAMPLE
		python SortTimePlot.py Results/Sorting/Hopfield_10_.txt
	
	HISTORY
		2016/01/28	Started
	"""
	me = "SortTimePlot.main: "
	t0 = time.time()
	
	try:
		argv[1]
	except IndexError:
		print main.__doc__
		raise IndexError(me+"Please read docstring")
	
	if os.path.isfile(argv[1]):
		timeplot(argv[1])
	elif os.path.isdir(argv[1]):
		plotall(argv[1])
	else:
		raise IOError("Input not a file or directory.")
	
	print me+"Execution",round(time.time()-t0,2),"seconds."
	return


##=============================================================================
def timeplot(datafile):
	"""
	"""
	me = "SortTimePlot.timeplot: "
	
	plotfile = datafile[:-4]+".png"
	Delta = findpars(datafile)
	
	npts = 500
	data = np.loadtxt(datafile, skiprows=4, unpack=True)
	data = data[:, ::int(data.shape[1]/npts)]
	
	t, ent, work, tran = data
	Dent, SSidx = flatline(ent, Delta)
	
	##-------------------------------------------------------------------------
	## Find average work rate and final entropy value
	## Assuming entropy is flat and work is linear
	
	S_fin = np.mean(ent[SSidx:])
	Wdot_evo = np.mean(work[:SSidx-int(npts/20)])/t[SSidx]
	Wdot_SS = np.mean(work[SSidx:]-work[SSidx])/(t[-1]-t[SSidx])
	annotext = "$\Delta S = %.2e$ \n $\dot W_{evo} = %.2e$ \n $\dot W_{SS} = %.2e$ \n $t_{SS} =  %.2e$"\
		% (S_fin, Wdot_evo, Wdot_SS, t[SSidx])
	
	##-------------------------------------------------------------------------
	## Plotting
	
	plt.clf()
	
	plt.plot(t, ent, "-", label="$S$")
	plt.plot(t, work, "-", label="$W$")

	plt.plot(t, 100*Dent, "--", label="$100\dot S$")
	plt.plot(t, 100*gaussian_filter1d(work,len(work)/100,order=1), "--", label="$100\dot W$")
	
	plt.vlines([t[SSidx-int(npts/20)],t[SSidx]],-4000, 500)
	plt.axvspan(t[0],t[SSidx-int(npts/20)], color="y",alpha=0.05)
	plt.axvspan(t[SSidx],t[-1], color="g",alpha=0.05)
	
	plt.ylim([-4000,500])
	plt.xlabel("$t$")
	plt.title("$\Delta=$"+str(Delta))
	plt.legend()
	plt.annotate(annotext,xy=(0.15,0.2),xycoords="figure fraction",fontsize=16)
	plt.grid()
	
	plt.savefig(plotfile)
	print me+"Plot saved to",plotfile
	# plt.show()
	
	return
	
##=============================================================================

def findpars(filename):
	start = filename.find("_") + 1
	Delta = float(filename[start:filename.find("_",start)])
	return Delta

##=============================================================================

def flatline(y,D):
	"""
	Iterates though array y and returns the index of first instance cloe to zero
	"""
	fac = 2 if D <=2.5 else 5
	y_conv = gaussian_filter1d(y,len(y)/200,order=1)
	for id, el in enumerate(y_conv):
		if abs(el)<2e-1 and id>len(y_conv)/fac:
			break
	return (y_conv, id)


##=============================================================================

def plotall():
	for datafile in glob.glob(argv[1]+"/*.txt"):
		timeplot(datafile)
	return

##=============================================================================
##=============================================================================
if __name__=="__main__":
	main()
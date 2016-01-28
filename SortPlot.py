
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d
from matplotlib import pyplot as plt
from sys import argv
import os, glob, time

##=============================================================================
def main():
	"""
	NAME
		SortPlot.py
	
	PURPOSE
	
	EXECUTION
	
	OUTPUT
	
	EXAMPLE
	
	HISTORY
	"""
	me = "SortPlot.main: "
	t0 = time.time()
	
	if os.path.isfile(argv[1]):
		timeplot(argv[1])
	elif os.path.isdir(argv[1]):
		dirplot(argv[1])
	else:
		raise IOError("Input not a file or directory.")
	
	print me+"Execution",round(time.time()-t0,2),"seconds."
	return


##=============================================================================
def timeplot(datafile):
	"""
	"""
	me = "SortPlot.timeplot: "
	
	plotfile = datafile[:-4]+".png"
	Delta = findpars(datafile)
	
	npts = 500
	data = np.loadtxt(datafile, skiprows=4, unpack=True)
	data = data[:, ::int(data.shape[1]/npts)]
	
	t, ent, work, tran = data
	
	plt.plot(t, ent, "-", label="$S$")
	plt.plot(t, work, "-", label="$W$")

	SSidx = flatline(ent)
	print SSidx, t[SSidx]
	plt.axvspan(t[0], t[SSidx-int(npts/10)],color="y",alpha=0.2)
	plt.axvspan(t[SSidx+int(npts/10)],t[-1],color="g",alpha=0.2)
	
	plt.ylim([-4000,500])
	plt.xlabel("$t$")
	plt.title("$\Delta=$"+str(Delta))
	plt.legend()
	plt.grid()
	
	plt.show()
	# plt.savefig(plotfile)
	print me+"Plot saved to",plotfile
	
	return
	
##=============================================================================

def findpars(filename):
	start = filename.find("_") + 1
	Delta = float(filename[start:filename.find("_",start)])
	return Delta

##=============================================================================

def flatline(y):
	"""
	Iterates though array y and returns the index of first instance cloe to zero
	"""
	y_conv = gaussian_filter1d(y,10,order=1)
	print y.shape, y_conv.shape
	print y_conv<-1e-5; exit()
	for id, el in enumerate(y_conv):
		if el<-1e-5:
			return id
	return None


##=============================================================================

def plotall(dir):
	for datafile in glob.glob(dir+"/*.txt"):
		timeplot(datafile)
	return

##=============================================================================
##=============================================================================
if __name__=="__main__":
	main()
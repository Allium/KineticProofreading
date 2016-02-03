
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
	data = np.loadtxt(datafile, skiprows=10, unpack=True)
	data = data[:, ::int(data.shape[1]/npts)]
	
	t, ent, work, trans_C, trans_I = data[[0,5,6,8,9]]
	N = int(data[[1,2,3,4],0].sum())
	del data
	
	ent /= N*np.log(2)	
	Dent, SSidx = flatline(ent, N, Delta)
	work /= np.abs(work[-1])
	
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
	
	plt.plot(t, ent, "b-", label="$S / N\ln2$")
	plt.plot(t, work, "g-", label="$W$")

	plt.axhline(y=SSS_theo(Delta),color="b",linestyle=":",linewidth=2, label="Hopfield")
	# plt.hlines(SSS_theo(Delta),t[0],t[-1],color="b",ls=":",linewidth=4)
	
	# mfac = 1
	# plt.plot(t, mfac*Dent, "b--", label="$%.0f\dot S$"%mfac)
	# plt.plot(t, mfac*gaussian_filter1d(work,len(work)/100,order=1), "g--", label="$%.0f\dot W$"%mfac)
	
	plt.vlines([t[SSidx-int(npts/20)],t[SSidx]],-2,1)
	plt.axvspan(t[0],t[SSidx-int(npts/20)], color="y",alpha=0.05)
	plt.axvspan(t[SSidx],t[-1], color="g",alpha=0.05)
	
	plt.ylim([-1.1,0.1])
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

def flatline(y,N,D):
	"""
	Iterates though array y and returns the index of first instance cloe to zero
	"""
	fac = 20#2 if D <=2.5 else 5
	y_conv = gaussian_filter1d(y,len(y)/200,order=1)
	for id, el in enumerate(y_conv):
		if abs(el)<0.2/N and id>len(y_conv)/fac:
			break
	return (y_conv, id)

##=============================================================================

def SSS_theo(D):
	"""
	Prediction for the SS entropy.
	Assumes a single Delta. See notes 24/01/2015.
	Normalised to be -1 when complete segregation (D high)
	"""
	DS = ( np.log(1+D*D) - D*D/(1+D*D)*np.log(D*D) - np.log(2) ) / np.log(2)
	return DS

##=============================================================================

def plotall(dirpath):
	for datafile in glob.glob(dirpath+"/*.txt"):
		timeplot(datafile)
	return

##=============================================================================

def findpars(filename):
	start = filename.find("_") + 1
	Delta = float(filename[start:filename.find("_",start)])
	return Delta
	
def read_header(datafile):
	"""
	Read header info from datafile.
	"""
	head = []
	f = open(datafile,'r')
	for i,line in enumerate(fh):
		if i is 10: break
		head += [line]
	f.close()
	return head
	
def get_headinfo(datafile):
	"""
	Hard-coded function to extract initial number and rates from header string.
	Better to use a dict but I can't be bothered.
	"""
	head = read_header(datafile)
	k_labels = head[3]
	k_values = np.fromstring(head[4], dtype=float, sep='\t')
	kp_labels = head[6]
	kp_values = np.fromstring(head[7], dtype=float, sep='\t')
	return [k_values, kp_values]

##=============================================================================
##=============================================================================
if __name__=="__main__":
	main()
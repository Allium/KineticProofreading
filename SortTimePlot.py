
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
		python SortTimePlot.py [filename]
		
	OUTPUT
		png plot saved in file directory
	
	EXAMPLE
		python SortTimePlot.py Results/Sorting/Hopfield_10_.txt
	
	BUGS
		get_headinfo parsing
	
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
	Delta = get_pars(datafile)
	
	npts = 500
	data = np.loadtxt(datafile, skiprows=10, unpack=True)
	data = data[:, ::int(data.shape[1]/npts)]
	
	t, ent, work, trans_C, trans_I = data[[0,5,6,8,9]]
	N = int(data[[1,2,3,4],0].sum())
	del data
	
	ent /= N*np.log(2)
	Dent, Sssid = flatline(ent,  N, Delta)
		
	##-------------------------------------------------------------------------
	## Find average work rate and final entropy value
	## Assuming entropy is flat and work is linear
	
	S_fin = np.mean(ent[Sssid:])
	Wdot_evo = np.mean(work[:Wssid-int(npts/20)])/t[Wssid]
	Wdot_SS = np.mean(work[Wssid:]-work[Wssid])/(t[-1]-t[Wssid])
	annotext = "$\Delta S = %.2e$ \n $\dot W_{evo} = %.2e$ \n $\dot W_{SS} = %.2e$ \n $t_{SS} =  %.2e$"\
		% (S_fin, Wdot_evo, Wdot_SS, t[Wssid])
	
	##-------------------------------------------------------------------------
	## Plotting
	
	plt.clf()
	
	plt.plot(t, ent, "b-", label="$S / N\ln2$")
	plt.plot(t, -work/work[-1], "g-", label="$W / W(tmax)$")

	plt.axhline(y=SSS_theo(Delta),color="b",linestyle=":",linewidth=2, label="Hopfield")
	
	# mfac = 1
	# plt.plot(t, mfac*Dent, "b--", label="$%.0f\dot S$"%mfac)
	# plt.plot(t, mfac*gaussian_filter1d(work,len(work)/100,order=1), "g--", label="$%.0f\dot W$"%mfac)
	
	plt.vlines([t[Wssid-int(npts/20)],t[Wssid]],-2,1)
	plt.axvspan(t[0],t[Wssid-int(npts/20)], color="y",alpha=0.05)
	plt.axvspan(t[Wssid],t[-1], color="g",alpha=0.05)
	
	plt.xlim(left=0.0)
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

def flatline(y,N,D, order=1, kerfrac=200, eps=0.2):
	"""
	Iterates though array y and returns the index of first instance cloe to zero.
	Filter width is set as fraction of total length.
	Note no division by time step.
	"""
	fac = 20#2 if D <=2.5 else 5
	y_conv = gaussian_filter1d(y,len(y)/kerfrac,order=order)
	for id, el in enumerate(y_conv):
		if abs(el)<eps/N and id>len(y_conv)/fac:
			break
	return (y_conv, id)

##=============================================================================

def SSS_theo(D):
	"""
	Prediction for the SS entropy in equilibrium.
	For Hopfield D->D^2.
	Assumes a single Delta. See notes 24/01/2015.
	Normalised to be -1 when complete segregation (D high)
	"""
	return ( np.log(1+D) - D/(1+D)*np.log(D) - np.log(2) ) / np.log(2)

	
def SSW_theo(k,D,i):
	"""
	Prediction for the SS work RATE.
	i=0: Hopfield; i=1: notfield.
	Assumes a single value of D and everything else is the same between networks.
	Note A "wants" to be in box 2.
	Unnormalised at the moment.
	"""
	me = "SortTimePlot.SSW_theo: "
	A2_ss = D*D/(1+D*D) if i==0 else D/(1+D)
	try:
		C_ss = A2_ss * k["A1B1"]*k["B1C1"] /\
			( k["B1C1"]*(k["C1A2"]+k["C1A1"]*D) +\
			  k["B1A1"]*D*(k["C1B1"]+k["C1A2"]+k["C1A1"]*D) +\
			  A2_ss*k["A1B1"]*(k["B1C1"]+k["C1B1"]+k["C1A2"]+k["C1A1"]*D) )
		return C_ss*(2*k["C1A2"]+(1+D)*k["C1A1"])
	except TypeError:
		return A2_ss * (3+D) / ( 1+3*D+D*D+A1_ss*(3+D) )
	except KeyError:
		raise KeyError(me+"Check k-values in file header.\n"+k.tostring())
	
##=============================================================================

def plotall(dirpath):
	for datafile in glob.glob(dirpath+"/*.txt"):
		timeplot(datafile)
	return

##=============================================================================

def get_pars(filename):
	start = filename.find("_") + 1
	Delta = float(filename[start:filename.find("_",start)])
	return Delta
	
def read_header(datafile):
	"""
	Read header info from datafile.
	"""
	head = []
	f = open(datafile,'r')
	for i,line in enumerate(f):
		if i is 10: break
		head += [line]
	f.close()
	return head
	
def get_headinfo(datafile):
	"""
	Hard-coded function to extract initial number and rates from header string.
	There is a bug here in the parsing -- sometimes misses an entry.
	"""
	head = read_header(datafile)
	k_labels = [i for i in head[3].split("\t")] + [j for j in head[4].split("\t")]
	k_values = np.concatenate([np.fromstring(head[4],dtype=float,sep="\t"),
		np.fromstring(head[7],dtype=float,sep="\t")])
	k_dict = dict((label,k_values[i]) for i,label in enumerate(k_labels))
	return k_dict

##=============================================================================
##=============================================================================
if __name__=="__main__":
	main()

import numpy as np
from scipy.signal import fftconvolve
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
	k = get_headinfo(datafile)

	npts = 500
	data = np.loadtxt(datafile, skiprows=10, unpack=True)
	data = data[:, ::int(data.shape[1]/npts)]

	t, ent, work, trans_C, trans_I = data[[0,5,6,8,9]]
	N = int(data[[1,2,3,4],0].sum())
	del data


	ent /= N*np.log(2)
	Sssid = flatline(ent)

	##-------------------------------------------------------------------------
	## Find average work rate and final entropy value
	## Assuming entropy is flat and work is linear

	S_fin = np.mean(ent[Sssid:])
	Wdot_evo = np.mean(work[:Sssid-int(npts/20)])/t[Sssid]
	Wdot_SS = np.mean(work[Sssid:]-work[Sssid])/(t[-1]-t[Sssid])
	annotext = "$\Delta S = %.2e$ \n $\dot W_{evo} = %.2e$ \n $\dot W_{SS} = %.2e$ \n $t_{SS} =  %.2e$"\
		% (S_fin, Wdot_evo, Wdot_SS, t[Sssid])

	Wssid = Sssid

	##-------------------------------------------------------------------------
	## Plotting

	plt.clf()

	plt.plot(t, ent, "b-", label="$S / N\ln2$")
	plt.plot(t, -work/work[-1], "g-", label="$W / W(tmax)$")

	expnt = 2.0
	if os.path.basename(datafile)[:3] == "Not": expnt = 1.0
	plt.axhline(y=SSS_theo(Delta**expnt),color="b",linestyle=":",linewidth=2, label="Hopfield")

	# mfac = 1
	# plt.plot(t, mfac*Dent, "b--", label="$%.0f\dot S$"%mfac)
	# plt.plot(t, mfac*gaussian_filter1d(work,len(work)/100,order=1), "g--", label="$%.0f\dot W$"%mfac)

	plt.vlines(SSt_theo(k)*N,-2,1,color="r",linestyle=":",lw=2)

	plt.vlines([t[Wssid]],-2,1)
	plt.axvspan(t[0],t[Wssid], color="y",alpha=0.05)
	plt.axvspan(t[Wssid],t[-1], color="g",alpha=0.05)

	plt.xlim(left=0.0,right=t[-1])
	plt.ylim([-1.1,0.1])
	plt.xlabel("$t$")
	plt.title("$\Delta=$"+str(Delta))
	plt.legend()
	# plt.annotate(annotext,xy=(0.15,0.2),xycoords="figure fraction",fontsize=16)
	plt.grid()

	plt.savefig(plotfile)
	print me+"Plot saved to",plotfile
	# plt.show()

	return

##=============================================================================

def flatline(y):
	sslvl = np.mean(y[int(0.75*y.size):])
	win = y.size/20
	y_conv = fftconvolve(y,np.ones(win)/win,mode="same")
	for id, el in enumerate(y_conv-sslvl):
		if el<0.05*np.abs(sslvl-y[0]): break
	return id

##=============================================================================

def SSS_theo(D, nu=1.0):
	"""
	Prediction for the SS entropy in equilibrium.
	For Hopfield D->D^2.
	Assumes a single Delta. See notes 24/01/2015.
	Normalised to be -1 when complete segregation (D high)
	"""
	D = D**nu
	return ( np.log(1+D) - D/(1+D)*np.log(D) - np.log(2) ) / np.log(2)


def SSW_theo(D,k,nu=1.0):
	"""
	Prediction for the SS work RATE.
	i=0: Hopfield; i=1: notfield.
	Assumes a single value of D and everything else is the same between networks.
	Note A "wants" to be in box 2.
	Unnormalised at the moment.
	"""
	me = "SortTimePlot.SSW_theo: "
	# A2_ss = D*D/(1+D*D) if i==0 else D/(1+D)
	A2_ss = D**nu/(1+D**nu)
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

def SSt_theo(k):
	"""
	See notes 26/02/2016
	"""
	den1 = k["C1A2"]*k["B1C1"]*k["A1B1"]
	num1 = k["B1A1"]*(k["C1B1"]+k["C1A1"]+k["C1A2"])+k["B1C1"]*(k["C1A1"]+k["C1A2"])
	try:
		den2 = k["C1'A2'"]*k["B1'C1'"]*k["A1'B1'"]
		num2 = k["B1'A1'"]*(k["C1'B1'"]+k["C1'A1'"]+k["C1'A2'"])+k["B1'C1'"]*(k["C1'A1'"]+k["C1'A2'"])
	except KeyError:
		den2 = k["C1A2p"]*k["B1C1p"]*k["A1B1p"]
		num2 = k["B1A1p"]*(k["C1B1p"]+k["C1A1p"]+k["C1A2p"])+k["B1C1p"]*(k["C1A1p"]+k["C1A2p"])
	return min(num1/den1,num2/den2)

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
	k_labels = [label.rstrip() for label in head[3].split("\t")] + [label.rstrip() for label in head[6].split("\t")]
	k_values = np.concatenate([np.fromstring(head[4],dtype=float,sep="\t"),
		np.fromstring(head[7],dtype=float,sep="\t")])
	k_dict = dict((label,k_values[i]) for i,label in enumerate(k_labels))
	return k_dict

##=============================================================================
##=============================================================================
if __name__=="__main__":
	main()
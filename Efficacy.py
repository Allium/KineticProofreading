
import numpy as np
from scipy import ndimage
from matplotlib import pyplot as plt
import optparse
from sys import argv
import os, glob
from time import time as sysT
from HopfieldPlotter import plot_acco
from Plotter import read_header, get_headinfo

def main():
	"""
	NAME:
		Efficacy.py	
	
	PURPOSE:
		Efficacy metrics for a Hopfield and control simulations, as a function of
			the discrimination factor Delta.
	
	EXECUTION:
		python Efficacy.py Directory PlotType Flags
	
	ARGUMENTS:
		Directory	Containes files to be analysed
		PlotType	Select which plot is to be output
	
	OPTIONS:
		PlotType
			2	Ratio of product rates
			3	Energy cost
			5	Time delay
	
	FLAGS
		-p --npoints	200		How many points to keep for plotting
		-s --showfig	False	
	
	OUTPUTS:
		Desired plots saved to data directory in .png format
	
	COMMENTS
	
	BUGS AND TODO:
		-- Use higher-order centered derivatives
	
	HISTORY:
		03/11/2015 Started CS
	"""
	me = "Efficacy.main: "
	t0 = sysT()
	
	## CLA check
	try:
		1/int(os.path.isdir(argv[1]))
		int(argv[2])
	except:
		print me+" input issue. Abort.\n\n",main.__doc__; exit()
	
	## Networks considered
	network = ["no-loop ","Hopfield "]
	
	## Options
	parser = optparse.OptionParser()	
	parser.add_option('-p','--npoints',
		dest="npoints", default=200, type="int")
	parser.add_option('-s','--show',
		dest="showfig", default=False, action="store_true")
	parser.add_option('-v','--verbose',
		dest="verbose", default=False, action="store_true")
	opt = parser.parse_args()[0]
	npoints = opt.npoints
	showfig = opt.showfig
	verbose = opt.verbose
	
	##-------------------------------------------------------------------------
	## DATA COLLECTION
	
	## Files in directory
	datafiles = np.sort(glob.glob("Results/*.txt"))
	numfiles = len(datafiles); assert numfiles%2==0
	
	Delta = np.zeros(numfiles/2)
	if int(argv[2])==2:
		ncurves = 2
		Cordi = np.zeros((numfiles/2,ncurves))
		Hordi = np.zeros((numfiles/2,ncurves))
	else:
		pass
	
	## Loop over files
	for i in range(numfiles/2):
		
		## Filenames
		Hfile = datafiles[i]
		Cfile = datafiles[i+numfiles/2]
		
		## Unpack data
		time, Cdata, Hdata, Ck, Hk = unpack_data(Cfile, Hfile, npoints, verbose)
		
		## Extract Delta
		Delta[i] = Ck[6]/Ck[1]
		assert Delta[i]==Hk[7]/Hk[4]==Hk[6]/Hk[1]
		
		## Values for ordinate
		if int(argv[2])==2:
			Cordi[i] = prod_rate(Cdata[[1,2]],ncurves,True)
			Hordi[i] = prod_rate(Hdata[[1,2]],ncurves,True)
		else:
			print me+"functionality not written yet."; exit()
	
	## Sort data
	sortind = np.argsort(Delta)
	Delta = Delta[sortind]
	Cordi = Cordi[sortind]
	Hordi = Hordi[sortind]
	
	##-------------------------------------------------------------------------
	## PLOTTING
	
	times = np.linspace(0,1,ncurves+2).round(2).astype("string")
	colors = ["k","b","r","g","y"]
	
	ax = plt.figure().add_subplot(111)
	for i in range(ncurves):
		ax.plot(Delta,Cordi[:,i], colors[i]+"--",label=network[0]+times[i+1]+"*tmax")
		ax.plot(Delta,Cordi[:,i], colors[i]+"x")
		ax.plot(Delta,Hordi[:,i], colors[i]+"-" ,label=network[1]+times[i+1]+"*tmax")
		ax.plot(Delta,Hordi[:,i], colors[i]+"o")
	plot_acco(ax, title="Ratio of product rates for "+network[0]+"and Hopfield networks",\
				xlabel="$\Delta$", ylabel="Incorrect / Correct Product Rates")	
	ax.set_xlim([Delta.min(),Delta.max()])
	## Save plot
	figfile = argv[1]+"/0Efficacy_ProdRateRatio.png"
	plt.savefig(figfile); print me+"Figure saved to",figfile
	
	## Wrap up
	if verbose: print me+"total execution time:",round(sysT()-t0,2),"seconds"
	if showfig: plt.show()
	
	return
	
##=============================================================================
def unpack_data(Cdatafile, Hdatafile, npoints, verbose):
	"""
	Read and pre-process data
	"""
	me = "Efficacy.unpack_data: "
	t0_unpack_data = sysT()
	
	## Read data from file
	try:
		Chead = read_header(Cdatafile)
		CN0, Cpars = get_headinfo(Chead)
		Cdata = np.loadtxt(Cdatafile,skiprows=8,unpack=True)
		Hhead = read_header(Hdatafile)
		HN0, Hpars = get_headinfo(Hhead)
		Hdata = np.loadtxt(Hdatafile,skiprows=8,unpack=True)
	except TypeError:
		print me+"file issue. Abort.\n\n"+main.__doc__
		exit()
	
	## Ensure parameters match between the two files
	assert Cdata.shape==Hdata.shape
	assert all(Cdata[0]==Hdata[0])
	assert CN0==HN0
		
	## Select relevant columns
	Cdata = Cdata[[0,4,8,9,10,11,12,13,14]]
	Hdata = Hdata[[0,4,8,9,10,11,12,13,14]]
	
	## Prune in time
	if Cdata.shape[1] > npoints:
		Cdata = Cdata[:,np.arange(0,Cdata.shape[1],Cdata.shape[1]/npoints)]
		Hdata = Hdata[:,np.arange(0,Hdata.shape[1],Hdata.shape[1]/npoints)]
		
	## Normalise
	Cdata[1:] /= CN0
	Hdata[1:] /= HN0
	
	## Define time array for clarity
	time = Cdata[0]
	
	if verbose: print me+"data preparation overhead",round(sysT()-t0_unpack_data,2),"seconds."
	
	return time, Cdata, Hdata, Cpars[1], Hpars[1]

##=============================================================================
def prod_rate(prod, ntimes, ratio=False):
	"""
	Calculate the ratio of product formation rates (incorrect / correct rate)
	over ntimes intervals.
	prod is an array whose last dimension hols the time series data.
	"""
	me = "Efficacy.prod_rate_ratio: "
	## Cut off trailing points
	prod = prod[:,:-(prod.shape[-1]%ntimes)]
	## Split into equal chunks of time
	chunks = np.array_split(prod, ntimes, axis=1)
	## Centered derivative
	prodrate = ndimage.convolve1d(chunks,[-1,0,1],axis=-1).mean(axis=-1)
	if ratio:	return prodrate[1]/prodrate[0]
	else:		return prodrate

##=============================================================================
if __name__ =="__main__":
	main()
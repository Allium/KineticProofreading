
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
		-n --ncurves	1		How many chunks of time to plot separately
		-s --showfig	False	
		-f --fit		True	Exponential fit
		-v --verbose	False	Print information to screen
	
	OUTPUTS:
		Desired plots saved to data directory in .png format
	
	COMMENTS
	
	BUGS AND TODO:
		-- The linear fit for energy rates doesn't work
	
	HISTORY:
		03/11/2015	Started CS
		05/11/2015	Product rate comparison plot
		11/11/2015	Linear fit for averaging rates
		14/11/2015	Seek / save data file
					Fit rate ratio vs Delta to exponential
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
	parser.add_option('-n','--ncurves',
		dest="ncurves", default=1, type="int")
	parser.add_option('-s','--show',
		dest="showfig", default=False, action="store_true")
	parser.add_option('--nofit',
		dest="nofit", default=False, action="store_true")
	parser.add_option('-v','--verbose',
		dest="verbose", default=False, action="store_true")
	opt = parser.parse_args()[0]
	npoints = opt.npoints
	ncurves = opt.ncurves
	showfig = opt.showfig
	fit		= not opt.nofit
	verbose = opt.verbose
	
	##-------------------------------------------------------------------------
	## DATA COLLECTION
	## Note C means control and H means Hopfield
	
	## Plot strings
	## Product rate ratios
	if int(argv[2])==2:
		plotit = "Ratio of product rates for "+network[0]+"and Hopfield networks"
		ylabel = "Incorrect / Correct Product Rates"
		figfile = argv[1]+"/0Efficacy_ProdRateRatio_"+str(ncurves)+".png"
	## Energy costs
	elif int(argv[2])==3:
		plotit = "Energy dissipation rate for "+network[0]+"and Hopfield networks"
		ylabel = "Energy dissipation rate"
		figfile = argv[1]+"/0Efficacy_EnergyRateRatio_"+str(ncurves)+".png"
	else:
		print me+"functionality not written yet. Abort."; exit()
	
	outfile = os.path.splitext(figfile)[0]+".npy"
	if verbose: print me+"plotting:",plotit,"at",ncurves,"times"
	
	## Look for existing file; if none exists, make calculations and save
	##---------------------------------
	## Try to find file
	try:
		data = np.load(outfile).T
		Delta = data[0][:,np.newaxis]
		Cordi = data[1:ncurves+1].T
		Hordi = data[ncurves+1:].T
		del data
		if verbose: print me+"precomputed datafile found",outfile
		
	except IndexError:
		if verbose: print me+"something wrong with the datafile. Generating afresh."
		raise IOError
		
	##---------------------------------
	## If not found, read datafiles in directory and construct data arrays
	except IOError:
		## Find files in directory
		datafiles = np.sort(glob.glob(argv[1]+"/*.txt"))
		numfiles = len(datafiles); assert numfiles%2==0; numfiles = numfiles/2
		if verbose: print me+"found",numfiles,"file pairs"
		
		## Initialise arrays
		Delta = np.zeros(numfiles)
		Cordi = np.zeros((numfiles,ncurves))
		Hordi = np.zeros((numfiles,ncurves))
		
		## Loop over files
		for i in range(numfiles):
			
			## Filenames
			Hfile = datafiles[i]
			Cfile = datafiles[i+numfiles]
			
			## Unpack data (my function below)
			## Time-series and k-rates
			Cdata, Hdata, Ck, Hk = unpack_data(Cfile, Hfile, npoints, verbose)
			
			## Extract Delta from rates in header; ensure that files match
			Delta[i] = Ck[6]/Ck[1]
			assert Delta[i]==Hk[7]/Hk[4]==Hk[6]/Hk[1]
			
			## Values for ordinate (using my function below)
			## Product rate
			if int(argv[2])==2:
				Cordi[i] = prod_rate(Cdata[[1,2]],ncurves,True)
				Hordi[i] = prod_rate(Hdata[[1,2]],ncurves,True)
			## Energy rate
			elif int(argv[2])==3:
				CErate = prod_rate(Cdata[3:].sum(axis=0),ncurves,False)
				HErate = prod_rate(Hdata[3:].sum(axis=0),ncurves,False)
				Cordi[i] = CErate
				Hordi[i] = HErate
			else:
				print me+"functionality not written yet. Abort."; exit()
		
		## Sort data in order of increasing Delta
		sortind = np.argsort(Delta)
		Delta = Delta[sortind][:,np.newaxis]; Cordi = Cordi[sortind]; Hordi = Hordi[sortind]
		if verbose: print me+"Deltas found:",Delta.flatten()
		
		## Save data
		np.save(outfile,np.hstack([Delta,Cordi,Hordi]))
		if verbose: print me+"data saved to",outfile
	## End except
	
	##-------------------------------------------------------------------------
	## PLOTTING
	
	times = np.linspace(0,1,2*ncurves+1)[1::2].round(2).astype("string")
	colors = ["k","b","r","g"]
	
	ax = plt.figure().add_subplot(111)
	## Loop over times at which evaluation happens;
	## Plot points and lines for control and Hopfield
	for i in range(ncurves):
		ax.plot(Delta,Cordi[:,i], colors[i]+"--",label=network[0]+times[i]+"*tmax")
		ax.plot(Delta,Cordi[:,i], colors[i]+"x")
		ax.plot(Delta,Hordi[:,i], colors[i]+"-" ,label=network[1]+times[i]+"*tmax")
		ax.plot(Delta,Hordi[:,i], colors[i]+"o")
		## Also plot exponential fits
		if int(argv[2])==2 and fit:
			fitX,fitC,mC = exp_fit(Delta,Cordi[:,i])
			fitX,fitH,mH = exp_fit(Delta,Hordi[:,i])
			ax.plot(fitX, fitC , "m:", linewidth=2, label="$\exp["+str(round(mC,2))+"\Delta]$")
			ax.plot(fitX, fitH , "m:", linewidth=2, label="$\exp["+str(round(mH,2))+"\Delta]$")
	plot_acco(ax, title=plotit, xlabel="$\Delta$", ylabel=ylabel)	
	ax.set_xlim([Delta.min(),Delta.max()])
	ax.set_ylim(bottom=0.0)
	if np.hstack([Cordi,Hordi]).flatten().max() > 2.0: ax.set_ylim(top=2.0)
	## Save plot
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
	## 0:time; 4:prod; 8:prod'
	Cdata = Cdata[[0,4,8,9,10,11,12,13,14]]
	Hdata = Hdata[[0,4,8,9,10,11,12,13,14]]
	
	## Prune in time
	if Cdata.shape[1] > npoints:
		Cdata = Cdata[:,np.arange(0,Cdata.shape[1],Cdata.shape[1]/npoints)]
		Hdata = Hdata[:,np.arange(0,Hdata.shape[1],Hdata.shape[1]/npoints)]
		
	## Normalise per particle
	Cdata[1:] /= CN0
	Hdata[1:] /= HN0
		
	if verbose: print me+"data preparation overhead",round(sysT()-t0_unpack_data,2),"seconds."
	
	return Cdata, Hdata, Cpars[1], Hpars[1]

##=============================================================================
def prod_rate(prod, ntimes, ratio=False):
	"""
	Calculate the ratio of product formation rates (incorrect / correct rate)
	over ntimes intervals.
	prod is an array whose last dimension holds the time series data.
	If we consider product rate ratio, first dimension is product numbers for
	correct and incorrect substrates. Second dimension is time series.
	If we consider energy, the first dimension is the time series.
	"""
	me = "Efficacy.prod_rate: "
	## Are we comparing incorrect and correct? If so, w=2; if not, w=1
	w = len(prod.shape)
	## Cut trailing time series points from end
	if (prod.shape[-1]%ntimes) != 0:
		prod = prod[:,:-(prod.shape[-1]%ntimes)]
	## Split into equal chunks of time
	chunks = np.array_split(prod, ntimes, axis=-1)
	## Derivative
	# prodrate = ndimage.convolve1d(chunks,[1,-1],axis=-1)
	prodrate = ndimage.filters.gaussian_filter1d(chunks,10,axis=-1,order=1)
	## Mean
	if w==1:	## Crude!
		prodrate = prodrate.mean(axis=-1).T
	## Fit linear
	else:
		times = np.array_split(np.arange(0,prod.shape[-1],1.0),ntimes)
		prodrate = np.array([np.polyfit(times[j],chunks[j][i],1)[0]\
			for i in range(w) for j in range(ntimes)]).reshape([w,ntimes])
	if ratio:	return prodrate[1]/prodrate[0]
	else:		return prodrate
	
##=============================================================================
def exp_fit(x,y):
	"""
	Make an exponential fit to points  y(x).
	Returns new x and y coordinates.
	Intercept should be 1.0 exactly; I allow it to vary here.
	"""
	fit = np.polyfit(x.flatten(),np.log(y.flatten()),1)
	fit_fn = np.poly1d(fit)
	X = np.linspace(0,x[-1],5*x.size)
	return X, np.exp(fit_fn(X)), fit[0]

##=============================================================================
if __name__ =="__main__":
	main()
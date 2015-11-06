
import numpy as np
import scipy as sp
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
from sys import argv
import optparse
from time import time as sysT
import os, warnings
from Plotter import read_header, get_headinfo

np.seterr(divide='ignore',invalid='ignore')
warnings.filterwarnings("ignore",
	"No labelled objects found. Use label='...' kwarg on individual plots.",
	UserWarning)

def main():
	"""
	NAME:
		HopfieldPlotter.py	
	
	PURPOSE:
		Generates a variety of plots comparing the kinetics of the Hopfield
			network to those of a "control" four-state network.
	
	EXECUTION:
		python HopfieldPlotter.py datafile1 datafile2 plot_type flags
	
	ARGUMENTS:
		file1	Time series of reactant concentrations for "control" network
		file2	Time series of reactant concentrations for Hopfield network
		plot_type	Select which plot is to be output
	
	OPTIONS:
		plot_type
			2	Full kinetics for both networks
			3	Product kinetics for both networks on one plot
			5	Ratio of products for both networks
			7	Ratio of product rates for both networks
			11	Time delay for the networks
			13	Energy consumed by the networks
			n	A combination of the above based on the prime factors of n
	
	FLAGS
		-p --npoints	200		How many points to keep for plotting
		-s --showfig	False	
	
	OUTPUTS:
		outfile#	Desired plot(s) saved to data directory in .png format
	
	COMMENTS
	
	BUGS AND TODO:
		-- Choose npoints more cleverly
		-- bounds_error in interpolation
		-- energy plot (13) -- not sure what is interesting
	
	HISTORY:
		12/10/2015 Started
		17/10/2015 Usability and axis labelling improved
		23/10/2015 Plotting time delays
		02/10/2015 Energy consumption multiplot
	"""
	me = "HopfieldPlotter.py: "
	t0 = sysT()
	
	##-------------------------------------------------------------------------
	## INPUT AND DATA READ
	
	## Read user input: data files and plot_type
	try:
		Cdatafile = argv[1]
		Hdatafile = argv[2]
		argv[3]
	except IndexError:
		print me+"input issue. Abort.\n\n"+main.__doc__
		exit()
	## Options
	parser = optparse.OptionParser()	
	parser.add_option('-p','--npoints',
		dest="npoints", default=200, type="int")
	parser.add_option('-s','--show',
		dest="showfig", default=False, action="store_true")
	opt = parser.parse_args()[0]
	showfig = opt.showfig
	
	## Output filename base
	plotfile = os.path.splitext(Hdatafile)[0]
	
	## Networks considered
	network = ["no-loop ","Hopfield "]
	
	## Unpack data from file
	try:
		Ehead = read_header(Cdatafile)
		EN0, Epars = get_headinfo(Ehead)
		Cdata = np.loadtxt(Cdatafile,skiprows=8,unpack=True)
		Hhead = read_header(Hdatafile)
		HN0, Hpars = get_headinfo(Hhead)
		Hdata = np.loadtxt(Hdatafile,skiprows=8,unpack=True)
	except TypeError:
		print me+"file issue. Abort.\n\n"+main.__doc__
		exit()
		
	##-------------------------------------------------------------------------
	## DATA PREPARARION

	## Prune for speed and memory
	if Cdata.shape[1] > opt.npoints:
		Cdata = Cdata[:,np.arange(0,Cdata.shape[1],Cdata.shape[1]/opt.npoints)]
		Hdata = Hdata[:,np.arange(0,Hdata.shape[1],Hdata.shape[1]/opt.npoints)]
		
	## Normalise
	Cdata[1:] /= EN0
	Hdata[1:] /= HN0
		
	## Ensure parameters match between the two files
	assert Cdata.shape==Hdata.shape
	assert all(Cdata[0]==Hdata[0])
	assert EN0==HN0
	print me+"simulation parameters\n",\
		network[0],"{k}\t",Epars[1],"\n",\
		network[1],"{k}\t",Hpars[1]
	
	## Define time array for clarity
	time = Cdata[0]
	
	##-------------------------------------------------------------------------
	## PLOTTING
	
	colours = ["k","b","g","r","c","m","y"]
	labels  = ["Correct ", "Incorrect "]
	
	## Plot concentrations separately for naive and for Hopfield
	if float(argv[3])%2==0:
		fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
		for i in range(1,5):
			ax1.plot(time,Cdata[i,:],  colours[i]+"-", label=labels[0]+str(i))
			ax1.plot(time,Cdata[i+4,:],colours[i]+"--",label=labels[1]+str(i))
		plot_acco(ax1, title=network[0],ylabel="Fraction $N/N_0$")	
		for i in range(1,5):
			ax2.plot(time,Hdata[i,:],  colours[i]+"-", label=labels[0]+str(i))
			ax2.plot(time,Hdata[i+4,:],colours[i]+"--",label=labels[1]+str(i))
		plot_acco(ax2, title=network[1])		
		## Save plot
		figfile = plotfile+"_2FullKinetics.png"
		plt.savefig(figfile); print me+"Figure saved to",figfile
		if showfig: plt.show()
	
	## Plot naive products against Hopfield products
	if float(argv[3])%3==0:
		fig = plt.figure(); ax = fig.add_subplot(111)
		ax.plot(time,Cdata[4,:], "k"+"-", label=network[0]+labels[0])
		ax.plot(time,Cdata[8,:], "k"+"--",label=network[0]+labels[1])
		ax.plot(time,Hdata[4,:], "b"+"-", label=network[1]+labels[0])
		ax.plot(time,Hdata[8,:], "b"+"--",label=network[1]+labels[1])
		plot_acco(ax, title="Product production: "+network[0]+"versus Hopfield")	
		## Save plot
		figfile = plotfile+"_3ProdKinetics.png"
		plt.savefig(figfile); print me+"Figure saved to",figfile
		if showfig: plt.show()
	
	## Plot ratio of products for naive and Hopfield product
	if float(argv[3])%5==0:
		ax = plt.figure().add_subplot(111)
		ax.plot(time,Cdata[8,:]/Cdata[4,:], "k"+"-",label=network[0])
		ax.plot(time,Hdata[8,:]/Hdata[4,:], "b"+"-",label=network[1])
		ylabel = "Incorrect / Correct Product Concentration"
		plot_acco(ax, title="Ratio of products for "+network[0]+"and Hopfield networks",\
					ylabel=ylabel)		
		## Save plot
		figfile = plotfile+"_5ProdRatio.png"
		plt.savefig(figfile); print me+"Figure saved to",figfile
		if showfig: plt.show()

	## Plot derivative of naive product against Hopfield product
	if float(argv[3])%7==0:
		ax = plt.figure().add_subplot(111)
		# ax.plot(time[:-1],np.diff(Cdata[8,:])/np.diff(Cdata[4,:]),"k"+"-",label=network[0])
		# ax.plot(time[:-1],np.diff(Hdata[8,:])/np.diff(Hdata[4,:]),"b"+"-",label=network[1])
		envelope_plot(time[:-1],np.diff(Cdata[8,:],n=1)/np.diff(Cdata[4,:],n=1),\
			ax=ax,winsize=5,line="k-",label=network[0])
		envelope_plot(time[:-1],np.diff(Hdata[8,:],n=1)/np.diff(Hdata[4,:],n=1),\
			ax=ax,winsize=5,line="b-",label=network[1])
		plt.legend(network, loc="upper left")
		# ax.plot(time,hopfield_prod(Hdata[1],Hpars[1][:-2])/hopfield_prod(Hdata[1],Hpars[1][[0,6,2,3,7,5]])\
			# ,"b--",label="prediction")
		ylabel = "Incorrect / Correct Product Formation Rate"
		plot_acco(ax, title="Rate of product formation for "+network[0]+"and Hopfield networks",\
					ylabel=ylabel)
		## Save plot
		figfile = plotfile+"_7ProdRateRatio.png"
		plt.savefig(figfile); print me+"Figure saved to",figfile
		if showfig: plt.show()
	
	## Plot delay in time for different schemes to achieve a concentration of product
	if float(argv[3])%11==0:
		ax = plt.figure().add_subplot(111)
		ax.plot(*time_delay(time,Cdata[4],Hdata[4]),color="b",label=None)
		ylabel = "Time delay, $\Delta t$"
		plot_acco(ax, title="Time delay for product formation between "+network[0]+"and Hopfield networks",\
					xlabel="Product concentration $[P]$",ylabel=ylabel)
		## Save plot
		figfile = plotfile+"_11DtvP.png"
		plt.savefig(figfile); print me+"Figure saved to",figfile
		if showfig: plt.show()
	
	## Plot energy usage: versus time and correct product, scheme difference, and all contributions
	## Want breakdown of energy steps
	if float(argv[3])%13==0:
		nrgC = Cdata[9:12,:]; nrgH = Hdata[9:12,:]
		sumC = nrgC.sum(axis=0); sumH = nrgH.sum(axis=0)
		fig,ax = plt.subplots(2,2)
		## Versus time
		ax[0,0].plot(time,sumC,"k-",label=network[0])
		ax[0,0].plot(time,sumH,"b-",label="Hopfield")
		plot_acco(ax[0,0], title="Total $E$", ylabel="Energy units per particle, $E/N_0$",\
			legloc="upper left")
		ax[1,0].plot(time,sumH-sumC)
		plot_acco(ax[1,0], title="Difference $\Delta E$", ylabel="Energy units per particle, $E/N_0$")
		## Versus correct product
		ax[0,1].plot(Cdata[4]/Cdata[4,-1],sumC,"k-",label=network[0])
		ax[0,1].plot(Hdata[4]/Hdata[4,-1],sumH,"b-",label="Hopfield")
		plot_acco(ax[0,1], title="Total $E$", xlabel="Correct product (fraction of final)",\
			ylabel="Energy units per particle, $E/N_0$", legloc="upper left")
		## All energy consumtion pathways
		ax[1,1].plot(time,nrgC[0],"r--",label=network[0]+"BC")
		# ax[1,1].plot(time,nrgC[1],"g--",label=network[0]+"CA")
		ax[1,1].plot(time,nrgC[2],"b--",label=network[0]+"CD")
		ax[1,1].plot(time,nrgH[0],"r-", label="Hopfield BC")
		ax[1,1].plot(time,nrgH[1],"g-", label="Hopfield CA")
		ax[1,1].plot(time,nrgH[2],"b-", label="Hopfield CD")
		plot_acco(ax[1,1], title="Individual steps' consumption", ylabel="Energy units per particle, $E/N_0$",\
			legloc="upper left")
		plt.tight_layout()
		fig.set_size_inches(10,10)
		## Save plots
		figfile = plotfile+"_13Energy.png"
		plt.savefig(figfile); print me+"Figure saved to",figfile
		if showfig: plt.show()
	
	## Wrap up
	if not showfig: print me+"Execution time",round(sysT()-t0,2),"seconds"
	
	return

	
##=============================================================================
def plot_acco(ax, **kwargs):
	"""
	Plot accoutrements.
	kwargs: title, subtitle, xlabel, ylabel, plotfile
	"""
	me = "HopfieldPlotter.plot_acco: "
	try: ax.set_title(kwargs["title"])
	except: pass
	try: ax.suptitle(kawrgs["subtitle"])
	except: pass
	try: ax.set_xlabel(kwargs["xlabel"], fontsize=14)
	except: ax.set_xlabel("Time, $t$", fontsize=14)
	try: ax.set_ylabel(kwargs["ylabel"], fontsize=14)
	except: ax.set_ylabel("Fraction $N/N_0$", fontsize=14)
	ax.grid(True)
	try: ax.legend(loc=kwargs["legloc"], fontsize=11)
	except: ax.legend(loc="best", fontsize=11)
	return
	

##=============================================================================
def hopfield_prod(S,rates,E0=1):
	"""
	Rate of producing product.
	Notation and assumptions follow those in Hopfield 1974.
	"""
	kp,k,mp,m,l,w = rates
	return (kp*mp*w*E0*S / ((k+mp)*(l+w)+kp*(mp+l+w)*S))

##=============================================================================
def time_delay(x,y1,y2):
	"""
	"""
	me = "HopfieldPlotter.time_delay: "	
	## Interpolate x as a function of y: x=g(y)
	g1 = sp.interpolate.interp1d(y1,x, bounds_error=False)
	g2 = interp1d(y2,x, bounds_error=False)
	## Array of y-ticks
	yarr = np.linspace(0,round(np.append(y1,y2).max()+0.05,1),200)
	## Delta x
	Dx = g2(yarr)-g1(yarr)	
	return yarr,Dx
	

##================================================
def envelope_plot(x, y, **kwargs):
	"""
	Plot smoothed mean with grey envelope.
	"""
	## Read kwargs
	try: ax = kwargs["ax"]
	except KeyError: ax = plt.gca()
	try: winsize = kwargs["winsize"]
	except KeyError: winsize = len(x)/10
	try: fill = kwargs["fill"]
	except KeyError: fill = "grey"
	try: line = kwargs["line"]
	except KeyError: line = "b-"
	
	# Coarsely chunk the data, discarding the last window if it's not evenly
	# divisible.
	numwin = x.size // winsize
	ywin = y[:winsize * numwin].reshape(-1, winsize)
	xwin = x[:winsize * numwin].reshape(-1, winsize)
	# Find the min, max, and mean within each window 
	ymin = ywin.min(axis=1); ymax = ywin.max(axis=1)
	ymean = ywin.mean(axis=1)
	xmean = xwin.mean(axis=1)

	fill_artist = ax.fill_between(xmean, ymin, ymax, color=fill, 
		                  edgecolor='none', alpha=0.5)
	line, = ax.plot(xmean, ymean, line)
	return fill_artist, line

##================================================
	
##=============================================================================
##=============================================================================
if __name__=="__main__":
	main()
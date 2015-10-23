
import numpy as np
from matplotlib import pyplot as plt
from sys import argv
from time import time as sysT
import os
from Plotter import read_header, get_headinfo

np.seterr(divide='ignore',invalid='ignore')

def main():
	"""
	NAME:
		HopfieldPlotter.py	
	
	PURPOSE:
		Generates a variety of plots comparing the kinetics of the Hopfield
			network to those of a "control" four-state network.
	
	EXECUTION:
		python HopfieldPlotter.py datafile1 datafile2 plot_type
	
	ARGUMENTS:
		file1	Time series of reactant concentrations for "control" network
		file2	Time series of reactant concentrations for Hopfiled network
		plot_type	Select which plot is to be output
	
	OPTIONS:
		plot_type
			2	Full kinetics for both networks
			3	Product kinetics for both networks on one plot
			5	Ratio of products for both networks
			7	Ratio of product rates for both networks
			n	A combination of the above based on the prime factors of n
	
	OUTPUTS:
		outfile#	Desired plot(s) saved to data directory in .png format
					If option integer is "exact", also shows plot on-screen
	
	BUGS, COMMENTS AND TODO:
		-- Choose npoints more cleverly
	
	HISTORY:
		12/10/2015 Started
		17/10/2015 Usability and axis labelling improved
	"""
	me = "HopfieldPlotter.py: "
	# from time import time
	t0 = sysT()
	
	## Read user input data file
	try:
		Edatafile = argv[1]
		Hdatafile = argv[2]
	except IndexError:
		print main.__doc__
		exit()
	
	## Output filename base
	plotfile = os.path.splitext(Hdatafile)[0]
	
	## Networks considered
	network = ["no-loop ","Hopfield "]
	
	## Unpack data from file
	try:
		Ehead = read_header(Edatafile)
		EN0, Epars = get_headinfo(Ehead)
		Edata = np.loadtxt(Edatafile,skiprows=8,unpack=True)
		Hhead = read_header(Hdatafile)
		HN0, Hpars = get_headinfo(Hhead)
		Hdata = np.loadtxt(Hdatafile,skiprows=8,unpack=True)
	except TypeError:
		print me+"input issue. Abort."
		exit()
	
	## Normalise
	Edata[1:] /= EN0
	Hdata[1:] /= HN0
		
	## Ensure parameters match between the two files
	assert Edata.shape==Hdata.shape
	assert all(Edata[0]==Hdata[0])
	assert EN0==HN0
	print me+"simulation parameters\n",\
		network[0],"{k}\t",Epars[1],"\n",\
		network[1],"{k}\t",Hpars[1]
	
	## Prune for speed
	npoints = 200
	if Edata.shape[1]>npoints:
		Edata = Edata[:,np.arange(0,Edata.shape[1],Edata.shape[1]/npoints)]
		Hdata = Hdata[:,np.arange(0,Hdata.shape[1],Hdata.shape[1]/npoints)]
	
	time = Edata[0]
	
##=============================================================================
	## PLOTTING
	
	colours = ["k","b","g","r","c","m","y"]
	labels  = ["Correct ", "Incorrect "]
	
	## Plot concentrations separately for naive and for Hopfield
	if float(argv[3])%2==0:
		fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
		for i in range(1,5):
			ax1.plot(time,Edata[i,:],  colours[i]+"-", label=labels[0]+str(i))
			ax1.plot(time,Edata[i+4,:],colours[i]+"--",label=labels[1]+str(i))
		plot_acco(ax1, title=network[0],ylabel="Fraction $N/N_0$")	
		for i in range(1,5):
			ax2.plot(time,Hdata[i,:],  colours[i]+"-", label=labels[0]+str(i))
			ax2.plot(time,Hdata[i+4,:],colours[i]+"--",label=labels[1]+str(i))
		plot_acco(ax2, title=network[1])		
		plt.savefig(plotfile+"_FullKinetics.png")
		print me+"Figure saved to",plotfile+"_FullKinetics.png"
	
	## Plot naive products against Hopfield products
	if float(argv[3])%3==0:
		fig = plt.figure(); ax = fig.add_subplot(111)
		ax.plot(time,Edata[4,:], "k"+"-", label=network[0]+labels[0])
		ax.plot(time,Edata[8,:], "k"+"--",label=network[0]+labels[1])
		ax.plot(time,Hdata[4,:], "b"+"-", label=network[1]+labels[0])
		ax.plot(time,Hdata[8,:], "b"+"--",label=network[1]+labels[1])
		plot_acco(ax, title="Product production: "+network[0]+"versus Hopfield")		
		plt.savefig(plotfile+"_ProdKinetics.png")
		print me+"Figure saved to",plotfile+"_ProdKinetics.png"
	
	## Plot ratio of products for naive and Hopfield product
	if float(argv[3])%5==0:
		ax = plt.figure().add_subplot(111)
		ax.plot(time,Edata[8,:]/Edata[4,:], "k"+"-",label=network[0])
		ax.plot(time,Hdata[8,:]/Hdata[4,:], "b"+"-",label=network[1])
		ylabel = "Incorrect / Correct Product Concentration"
		plot_acco(ax, title="Ratio of products for "+network[0]+"and Hopfield networks",\
					ylabel=ylabel)		
		plt.savefig(plotfile+"_ProdRatio.png")
		print me+"Figure saved to",plotfile+"_ProdRatio.png"

	## Plot derivative of naive product against Hopfield product
	if float(argv[3])%7==0:
		ax = plt.figure().add_subplot(111)
		ax.semilogy(time[:-1],np.diff(Edata[8,:])/np.diff(Edata[4,:]),"k"+"-",label=network[0])
		ax.semilogy(time[:-1],np.diff(Hdata[8,:])/np.diff(Hdata[4,:]),"b"+"-",label=network[1])
		ylabel = "Incorrect / Correct Product Formation Rate"
		plot_acco(ax, title="Rate of product formation for "+network[0]+"and Hopfield networks",\
					ylabel=ylabel)
		plt.savefig(plotfile+"_DtProdRatio.png")
		print me+"Figure saved to",plotfile+"_DtProdRatio.png"
	
	## Wrap up
	print me+"Execution time",round(sysT()-t0,2),"seconds"
	if any(np.array([2,3,5,7]) == float(argv[3])): plt.show()
	
	return

	
##=============================================================================
def plot_acco(ax, **kwargs):
	"""
	Plot accoutrements.
	kwargs: title, subtitle, ylabel, plotfile
	"""
	me = "HopfieldPlotter.plot_acco: "
	try: ax.set_title(kwargs["title"])
	except: pass
	try: ax.suptitle(kawrgs["subtitle"])
	except: pass
	ax.set_xlabel("Time, $t$", fontsize=14)
	try: ax.set_ylabel(kwargs["ylabel"], fontsize=14)
	except: ax.set_ylabel("Fraction $N/N_0$", fontsize=14)
	ax.grid(True)
	ax.legend(loc="center right", fontsize=11)			
	return
	

##=============================================================================
def hopfield_error():
	"""
	Rate of producing incorrect versus correct product.
	Notation and assumptions follow those in Hopfield 1974.
	-- same concentration of both substrates.
	-- discrimination in off-rates
	-- ...
	"""

##=============================================================================
##=============================================================================
if __name__=="__main__":
	main()
	
	
	

import numpy as np
from matplotlib import pyplot as plt
from sys import argv
import os

def main():
	"""
	NAME:
		plotter.py	
	
	PURPOSE:
		A simple plotter for reaction kinetics. Compares simulated population with exact solution.
	
	EXECUTION:
		python plotter.py <datafile> <arg>
	
	INPUTS:
		<datafile>	three arrays: time, simulated data, exact solution
	
	ARGUMENTS:
		argv[2]:	forces rows for predicted values
	
	OUTPUTS:
		<datafile>.png	plot
	
	BUGS, COMMENTS AND TODO:
		-- Other reading methods need update
		-- Horrible scraping...
	
	HISTORY:
		25/09/2015 Started
		02/10/2015 Arbitrary number of inputs
		12/10/2015 Scrape file header for relevant information
	"""
	
	## Read user input data file
	try: datafile = argv[1]
	except IndexError: print main.__doc__; exit()
	
	## Generate output filename
	plotfile = os.path.splitext(datafile)[0]+".png"
	
	## Inclusion of exact result is an option
	try:
		pre = argv[2]=="pre"
	except IndexError:
		pre = False
	
	## Unpack arrays from file: time, simulation, exact
	## Assuming column headers in first row
	try:
		data = np.load(datafile).T
	except:
		try:
			head = read_header(datafile)
			N0, plotstring = get_headinfo(head)
			data = np.loadtxt(datafile,skiprows=8,unpack=True)
		except TypeError:
			print "haegd"
			exit()
			data = np.genfromtxt(datafile,dtype="string",skip_header=8,unpack=True).astype("float")
	
	## Prune for speed
	npoints = 200
	if data.shape[1]>npoints:	data = data[:,np.arange(0,data.shape[1],data.shape[1]/npoints)]
	
	## Number of particles
	npar = (data.shape[0]-1)
	t = data[0]
	if pre:
		n_pre = data[1:npar+1]/N0
		n_sim = data[npar+1:]/N0
		npar /= 2
	else:
		n_sim = data[1:]/N0
	del data
			
	## Plot data
	colours = ["k","b","g","r","c","m","y"]
	for i in xrange(npar):
		if pre: plt.plot(t,n_pre[i],colours[i]+"-",label=str(i)+" prediction")
		plt.plot(t,n_sim[i],colours[i]+"x",markersize=10,label=str(i)+" simulation")
	
	## Plot accoutrements
	plt.ylim([0.0,1.0])
	plt.title("Population Dynamics")
	try: plt.suptitle(str(plotstring))
	except: pass
	plt.xlabel("Time, $t$"); plt.ylabel("Fraction remaining, $N/N_0$")
	plt.grid(True)
	plt.legend()
	
	## Output
	plt.savefig(plotfile)
	print "plotter.py: plot saved as ./"+plotfile
	plt.show()
	
	return

	
##=============================================================================
def read_header(datafile):
	"""
	Shitty function to read the header info from datafile.
	"""
	head = []
	fh = open(datafile,'r')
	for i,line in enumerate(fh):
		if i is 8: break
		head += [line]
	fh.close()
	return head
	
def get_headinfo(head):
	"""
	Hard-coded function to extract initial number and rates from header string.
	Better to use a dict but I can't be bothered.
	"""
	N_labels = head[1][:-1]
	N_values = np.fromstring(head[2][1:-2], dtype=float, sep=',')
	k_labels = head[3][:-1]
	k_values = np.fromstring(head[4][1:-2], dtype=float, sep=',')
	headinfo = N_values.sum(), [N_values,k_values]
	return headinfo
	
##=============================================================================
##=============================================================================
if __name__=="__main__":
	main()
	
	
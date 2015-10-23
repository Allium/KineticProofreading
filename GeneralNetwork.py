
import numpy as np
import time

def main(datafile=None):
	"""
	NAME:		GeneralNetwork.py	
	
	PURPOSE:	Simulate fully-connected reaction paths.
	
	EXECUTION:	python GeneralNetwork.py
	
	INPUTS:		None
	
	OUTPUTS:	datafile.txt	Time series of data (sensible density)
	
	BUGS, COMMENTS AND TODO:
	-- Express t in units of kmin
	
	HISTORY:	03/10/2015 Started
	"""
	
	## Initial populations of different species, total number, number of species
	Ni0 = np.array([1e4,1e4,1e3,5e2],dtype="int")
	N0  = Ni0.sum()	
	nspec = Ni0.size
	
	## Particles are represented as an array of length N0.
	## Each element is a particle. The value denotes the particle identity.
	parts = np.hstack(np.append([],[j*np.ones(Ni0[j]) for j in xrange(nspec)]))
	
	## Reaction rates
	## ks format: [0->1,0->2,...,0->n,1->2,...,1->0,...,n->0,...,n->n-1]
	ks = np.array([0.01,0.05,0.0,0.01,0.0,0.01,0.01,0.02,0.01,0.0,0.0,0.04])
	assert ks.size == nspec*(nspec-1)
	
	## Outfile name
	if datafile is None: datafile = "./Results/GeneralNetwork"+str(nspec)
	
	## t is an integer counter
	tmax = int(1/ks[ks!=0.0].min())
	## Array to store particle numbers at each t
	numparts = np.zeros([tmax,nspec+1])
	numparts[:,0] = np.arange(0,tmax)
	
	# np.random.seed(10)
	t0 = time.time()
	for t in xrange(tmax):
		## Populations
		ptloc = np.array([parts == j for j in xrange(nspec)])
		## Random array
		R = np.random.rand(N0)
		
		## Particle switching: probability of advancing [1,2,...,n-1] for each particle
		## Dimensions (nspec-1)*nparts
		pswch = np.array( [np.outer(ks[j*(nspec-1):(j+1)*(nspec-1)],ptloc[j]) for j in xrange(nspec)] ).sum(axis=0)
		
		## Change in particle identity
		Dpart = 1*(R<=pswch[0]) + np.array( [(j+1)*( (pswch[j-1]<=R)*(R<=pswch[j])) for j in xrange(1,nspec-1)] ).sum(axis=0)
		parts = (parts + Dpart) % nspec
		
		## Store values
		numparts[t,1:] = [ptloc[j].sum() for j in xrange(nspec)]
	
	print "ChemSim.py: Loop",round(time.time()-t0,1),"seconds for",N0,"particles and",tmax,"steps"

	## Output
	np.save(datafile,numparts/N0)
	print "ChemSim.py: Data saved to",datafile+".npy"

	return

##=======================================================================================
if __name__=="__main__":
	main()

import numpy as np
import time

def main():
	"""
	NAME:		LinearNetwork.py	
	
	PURPOSE:	Simulate linear (loop-free) reaction paths.
	
	EXECUTION:	python LinearNetwork.py
	
	INPUTS:		None
	
	OUTPUTS:	datafile.txt	Time series of data
	
	BUGS, COMMENTS AND TODO:
	-- Express t in units of kmin
	
	HISTORY:	02/10/2015 Started
	"""
	
	## Initial populations of different species, total number, number of species
	Ni0 = np.array([1e4,1e4,1e3,5e3],dtype="int")
	N0  = Ni0.sum()	
	nspec = Ni0.size
	
	## Particles are represented as an array of length N0.
	## Each element is a particle. The value denotes the particle identity.
	parts = np.hstack(np.append([],[j*np.ones(Ni0[j]) for j in xrange(nspec)]))
	
	## Reaction rates; LINEAR; non-zero -> grand loop
	## Particle 0 rate forward, rate back, Particle 1, ...
	kf = np.array([0.04,0.02,0.06,0.0])
	kb = np.array([0.0,0.02,0.03,0.01])
	assert kf.size == kb.size == nspec
	
	## Outfile name
	datafile = "./Results/Linear"+str(nspec)
	
	## t is an integer counter
	tmax = int(10/np.array([kf,kb]).max())
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
		
		## Particle switching: probability of forward [0] and back [1] for each particle
		pswch = np.array([[kf[j]*ptloc[j],kb[j]*ptloc[j]] for j in xrange(nspec)]).sum(axis=0)
		Dpart = (R<=pswch[0]) - (R>=1-pswch[1])
		parts += Dpart
		parts %= nspec
		
		## Store values
		numparts[t,1:] = [ptloc[j].sum() for j in xrange(nspec)]
	
	np.save(datafile,numparts/N0)
	print "ChemSim.py: Data saved to",datafile+".npy"
	print "ChemSim.py: Loop",round(time.time()-t0,1),"seconds for",N0,"particles and",tmax,"steps"

	return

##=======================================================================================
if __name__=="__main__":
	main()
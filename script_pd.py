import numpy as np
import os, shutil, glob, time
from sys import argv

me = "script_pd: "
t0 = time.time()

##=== Edit this section

A1B = np.arange(0.01,0.04,0.01)
BA1 = np.arange(0.003,0.02,0.004)
BC  = np.arange(0.003,0.04,0.004)
CA1 = np.arange(0.003,0.02,0.004)

##===

klist = [[a1b,ba1,bc,ca1] for a1b in A1B for ba1 in BA1 for bc in BC for ca1 in CA1]
print me+str(len(klist))+" points. Time estimate: ~"+str(round(len(klist)*25.0/3600,1))+" hours."

os.system("javac JHopfield.java")
for D in [10.0]:
	outdir = "./Results/SymChannel/Phase/D"+str(D)+"_2/"
	if not os.path.isdir(outdir):
		os.makedirs(outdir)
		print me+"Created directory "+outdir
	shutil.copy("./script_pd.py",outdir)
	print me+"Saved a copy of this script in "+outdir
	for i,ks in enumerate(klist):
		ks = " ".join([str(k) for k in ks])
		t1 = time.time()
		print me+"Delta =",D,"  i =",i+1,"of",len(klist)
		os.system("java JHopfield "+str(D)+" h c "+ks+" "+str(i+1000))
		print me+str([i+1]),round(time.time()-t1,1),"seconds"
print me+"Total execution time",round((time.time()-t0)/3600,1),"hours for",len(klist),"points."
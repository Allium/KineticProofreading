import numpy as np
import os, shutil, glob, time
from sys import argv

me = "script_pd: "

os.system("javac JHopfield.java")

A1B = np.arange(0.005,0.04,0.01)
BA1 = np.arange(0.001,0.02,0.004)
BC  = np.arange(0.001,0.04,0.004)
CA1 = np.arange(0.001,0.02,0.004)

klist = [[a1b,ba1,bc,ca1] for a1b in A1B for ba1 in BA1 for bc in BC for ca1 in CA1]

print me+"Points:",len(klist),"  Time: ~",len(klist)*30.0/60/60,"hours."

for D in [10.0,15.0]:
	outdir = "./Results/SymChannel/Phase/D"+str(D)+"/"
	if not os.path.isdir(outdir):
		os.makedirs(outdir)
		print me+"Created directory "+outdir
	shutil.copy("./script_pd.py",outdir)
	print me+"Saved a copy of this script in "+outdir
	for i,ks in enumerate(klist):
		ks = " ".join([str(k) for k in ks])
		t0 = time.time()
		print me+"Delta =",D,"  i =",i+1,"of",len(klist)
		os.system("java JHopfield "+str(D)+" h c "+ks+" "+str(i))
		print me,i+1,round(time.time()-t0,1),"seconds"
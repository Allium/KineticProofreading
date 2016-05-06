import numpy as np
import os, shutil, glob, time
from sys import argv

me = "script_pd: "

os.system("javac JHopfield.java")

A1B = [0.04,0.01]
BA1 = [0.02,0.01,0.005,0.001]
BC  = [0.02,0.01,0.005,0.001]
CA1 = [0.02,0.01,0.005,0.001]
# A1B = [0.4]
# BC = [0.0005]
# BA1 = [0.02]
# CA1 = [0.005]

klist = [[a1b,ba1,bc,ca1] for a1b in A1B for ba1 in BA1 for bc in BC for ca1 in CA1]

for D in [15.0]:
	outdir = "./Results/SymChannel/Phase/D"+str(D)+"/"
	if not os.path.isdir(outdir):
		os.makedirs(outdir)
		print me+"Created directory "+outdir
	shutil.copy("./script_pd.py",outdir)
	print me+"Saved copy of script in "+outdir
	exit()
	for i,ks in enumerate(klist):
		ks = " ".join([str(k) for k in ks])
		t0 = time.time()
		print me+"Delta =",D,"\ti =",i+1,"of",len(klist)
		os.system("java JHopfield "+str(D)+" h c "+ks+" "+str(i))
		print me,i+1,round(time.time()-t0,1),"seconds"
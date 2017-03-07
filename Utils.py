me0 = "Utils"
import numpy as np
import os, time, subprocess
from datetime import datetime

"""
NAME
Utils.py

PURPOSE
Supporting functions for kinetic proofreding project.

EXECUTION
None
"""


##================================================
## PLOTTING
##================================================


fs = {"fsa":30,"fsl":26,"fst":20,"fsn":26,"figsize":(10,10),"saveext":"pdf"}
# fs = [14,12,14]
# figsize = (4,4)

def set_mplrc(fs):
	"""
	Set MPL defaults
	"""
	
	import matplotlib as mpl
	from cycler import cycler
	
	## Number format
	mpl.rc("axes.formatter", limits=(-3,3))
	
	## Lines
	mpl.rc("lines", linewidth=2, markersize=8)
#	mpl.rc("axes", prop_cycle=cycler('color',["348ABD","7A68A6","A60628","467821","CF4457","188487","E24A33"])) ## Not working

	## Labels and legend
	mpl.rcParams["xtick.labelsize"] = fs["fsn"]
	mpl.rcParams["ytick.labelsize"] = fs["fsn"]
	mpl.rc("axes", labelsize=fs["fsa"])
	mpl.rc("legend", fontsize=fs["fsl"], fancybox=True)#framealpha=0.5, 
	
	## Font
	mpl.rcParams['font.family'] = 'serif'  
	mpl.rcParams['font.serif'] = ['Computer Modern Roman']  
	mpl.rcParams['text.usetex'] = True
	## DFM
	# rcParams["font.family"] = "sans-serif"
	# rcParams["font.sans-serif"] = ["Computer Modern Sans"]
	# rcParams["text.usetex"] = True
	# rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
	
	## Figure properties
	mpl.rc("figure", figsize=fs["figsize"])
	mpl.rc("savefig", format="jpg", dpi=200)
	
	return

	
##==========================================
## INPUT / OUTPUT
##==========================================

def filename_par(filename, searchstr):
	"""
	Scrape filename for parameter
	"""
	start = filename.find(searchstr) + len(searchstr)
	finish = start + 1
	while unicode(filename[start:].replace(".",""))[:finish-start].isnumeric():
		finish += 1
	return float(filename[start:finish-1])


## ====================================================================

def check_path(histfile, vb):
	"""
	Check whether directory exists; and if existing file will be overwritten.
	"""
	me = "Utils.check_path: "
	if os.path.isfile(histfile):
		raise IOError(me+"file",histfile,"already exists. Not overwriting.")
	try:
		assert os.path.isdir(os.path.dirname(histfile))
	except AssertionError:
		os.mkdir(os.path.dirname(histfile))
		if vb: print me+"Created directory",os.path.dirname(histfile)
	return
	
def create_readme(histfile, vb):
	"""
	If no readme exists, make one.
	NOTE commit is the LAST COMMIT -- maybe there have been changes since then.
	Assumes directory exists.
	"""
	me = "Utils.create_readme: "
	readmefile = os.path.dirname(histfile)+"/README.txt"
	try:
		assert os.path.isfile(readmefile)
	except AssertionError:
		now = str(datetime.now().strftime("%Y-%m-%d %H.%M"))
		commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'])
		header = "Time:\t"+now+"\nCommit hash:\t"+commit+"\n\n"
		with open(readmefile,"w") as f:
			f.write(header)
		if vb: print me+"Created readme file "+readmefile
	return

## ====================================================================

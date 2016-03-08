
import numpy as np
from matplotlib import pyplot as plt
import glob
from sys import argv
from SortTimePlot import get_pars, flatline

def main():
	dir_name = argv[1]
	errorplot(argv[1])
	plt.savefig(dir_name+"Error_Rates.png")
	plt.show()
	return

def errorplot(dir_name, ax=None, fits=[1.4,0.9]):

	file_list = glob.glob(dir_name+"/Hop*.txt")        
	file_listN = glob.glob(dir_name+"/Notfield*.txt")

	file_list = np.sort(file_list)
	file_listN = np.sort(file_listN)

	numfiles = len(file_list)
	numfilesN = len(file_listN)

	Delta = np.zeros(numfiles)
	DeltaN = np.zeros(numfilesN)

	Error = np.zeros(numfiles)
	Errorb = np.zeros(numfiles)
	for i in range(numfiles):
			Delta[i] = get_pars(file_list[i])
			data = np.loadtxt(file_list[i], skiprows=10).T
			N = data[[1,2,3,4],0].sum()
			idx = int(.25*data.shape[1])
			idxb = int(.75*data.shape[1])
			Error[i] = data[12,idx:].mean()
			Errorb[i] = data[12, idxb:].mean()


	ErrorN = np.zeros(numfilesN)
	ErrorNb = np.zeros(numfilesN)
	for j in range(numfilesN):
			DeltaN[j]= get_pars(file_listN[j])
			dataN = np.loadtxt(file_listN[j], skiprows=10).T
			Nn = dataN[[1,2,3,4],0].sum()
			idxN = int(.25*dataN.shape[1])
			idxbN = int(.75*dataN.shape[1])
			ErrorN[j]= dataN[12, idxN:].mean()
			ErrorNb[j] = dataN[12, idxbN:].mean()

	srtidx = np.argsort(Delta)
	srtidxb = np.argsort(Delta)

	srtidxN = np.argsort(DeltaN)
	srtidxbN = np.argsort(DeltaN)

	Delta = Delta[srtidx]
	Error = Error[srtidx]
	Errorb = Errorb[srtidxb]

	DeltaN = DeltaN[srtidxN]
	ErrorN = ErrorN[srtidxN]
	ErrorNb = ErrorNb[srtidxbN]
	
	if ax == None:
		ax = plt.gca()
		ax.set_title("I./C.")
		ax.set_xlabel("$\Delta$")
		ax.set_ylabel("Hopfield's Error Rate versus $\Delta$")
		fsl = 12
	else:
		ax.set_ylabel("Error Rate Ratio")
		fsl = 6
		
		

	#       plt.plot(Delta, Error, "rx-", markersize = 10, label = "0.25 Data Hopfield")
	ax.plot(Delta, Errorb, "bo", label = "Hopfield")
	ax.plot(Delta, Delta**(-fits[0]), "b--", label = "$\Delta^{-"+str(fits[0])+"}$")

	#        plt.plot(DeltaN, ErrorN, "bx-", markersize = 10, label = "0.25 Data Notfield")
	ax.plot(DeltaN, ErrorNb, "ro", label = "Notfield")
	ax.plot(Delta, Delta**(-fits[1]), "r--", label = "$\Delta^{-"+str(fits[1])+"}$")

	ax.plot(Delta, Delta**(-2), "m-", label = "Ideal $\Delta^{-2}$")

	ax.set_yscale('log')

	ax.grid(True)
	ax.legend(loc="lower left",fontsize=fsl)

	return

if __name__ == "__main__":
	main()

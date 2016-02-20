
import numpy as np
from matplotlib import pyplot as plt
import glob
from sys import argv
from SortTimePlot import get_pars, flatline

def main():
        dir_name = argv[1]
        errorplot(argv[1])
        return

def errorplot(dir_name):
        file_list = glob.glob(dir_name+"/Hop*.txt")        
        file_listN = glob.glob(dir_name+"/Not*.txt")
        file_list = np.sort(file_list)
        file_listN = np.sort(file_listN)

        plotfile = "Error_Rates.png"

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
                Nn = data[[1,2,3,4],0].sum()
                idxN = int(.25*dataN.shape[1])
                idxbN = int(.75*dataN.shape[1])
                ErrorN[j]= data[12, idxN:].mean()
                ErrorNb[j] = data[12, idxbN:].mean()

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
        
	plt.title("I./C.")
	plt.xlabel("$\Delta$")
	plt.ylabel("Hopfield's Error Rate versus $\Delta$")
	
        plt.plot(Delta, Error, "rx-", markersize = 10, label = "0.25 Data Hopfield")
        plt.plot(Delta, Errorb, "ro-", label = "0.75 Data Hopfield")

        plt.plot(DeltaN, ErrorN, "bx-", markersize = 10, label = "0.25 Data Notfield")
        plt.plot(DeltaN, ErrorNb, "bo-", label = "0.75 Data Notfield")

        print Error
        print Errorb
        print ErrorN
        print ErrorNb
        
        
 #       plt.plot(Delta, 1/(1+Delta*Delta), "b-", label = "Delta Squared")
        plt.grid(True)
        plt.legend()

       
        plt.savefig(plotfile)
        plt.show()
                   
        return

if __name__ == "__main__":
        main()

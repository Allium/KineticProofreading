
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
        
        file_list = glob.glob(dir_name+"/*.txt")
        file_list = np.sort(file_list)

        numfiles = len(file_list)
        Delta = np.zeros(numfiles)
        Error = np.zeros(numfiles)
        for i in range(numfiles):
                Delta[i] = get_pars(file_list[i])
                data = np.loadtxt(file_list[i], skiprows=10).T
                N = data[[1,2,3,4],0].sum()
                idx = int(0.75*data.shape[1])
                Error[i] = data[11,idx:].mean()

        srtidx = np.argsort(Delta)
        Delta = Delta[srtidx]
        Error = Error[srtidx]
        
	plt.title("Incorrect over Correct")
	plt.xlabel("Delta")
	plt.ylabel("Average Incorrect over Correct")
        plt.plot(Delta, Error)
        plt.grid(True)
#        plt.legend()

       
#        plt.savefig(plotfile)
        plt.show()
                   
        return

if __name__ == "__main__":
        main()

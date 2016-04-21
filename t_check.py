from SymPlot import SSt_theo

D = 5.0
k = {"A1B1": 0.01, "B1A1": 0.005, "C1A1": 0.01, "C1B1": 0.002}

print "Predicted SS time",round(SSt_theo(D, k),0), "times number of particles"
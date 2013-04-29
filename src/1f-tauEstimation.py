# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 15:47:52 2013

@author: Milad H. Mobarhan
"""


from pylab import *
import percolationLibrary as PL

close("all")

exponentList = range(4,10)
L = zeros(len(exponentList))

kIndex=0
for k in exponentList:
    L[kIndex] = 2**k 
    kIndex+=1

Lx = L
Ly = L


###If 1D:
#    pc=0.99
pc = 0.59275
nSamples = 10000
idString = "nsamples" + str(nSamples)

for LIndex in range(len(L)):
    sValues,nValues = PL.numberDensity(nSamples, Lx[LIndex], Ly[LIndex],pc)
        
    plot(log10(sValues[:-1]), log10(nValues),label = "L = "+ ("%.2f"% Lx[LIndex]))
#    savetxt("../results/1f/data-" + idString + ".dat", [sValues[:-1], nValues])
    print Lx[LIndex]
        
xlabel(r"$\log_{10}[s]$")
ylabel(r"$\log_{10}[n(s,p_c)]$")
legend(prop={'size':11})
xlim(xmin=1)
grid()
savefig("../results/1f/n-vs-s-" + idString + ".pdf")
           




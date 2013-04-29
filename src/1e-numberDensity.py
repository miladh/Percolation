# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 11:22:41 2013

@author: Milad H. Mobarhan
"""

from pylab import *
from scipy.ndimage import measurements
import percolationLibrary as PL

close("all")

Lx = 500
Ly = Lx
npValues = 5
nSamples = 1000

idString = "L" + str(Lx) +"x"+str(Ly) + "-nsamples" + str(nSamples)

for i in range(1):
    if i == 0:
        pValues  = linspace(0.5, 0.59, npValues)
        idString += "-frombelow"
    else:
        pValues  = linspace(0.60, 0.64, npValues)
        idString += "-fromabove"
        
    for pIndex in range(len(pValues)):
        p = pValues[pIndex]    
        
        sValues,nValues = PL.numberDensity(nSamples, Lx, Ly,p)
        
        plot(log10(sValues[:-1]), log10(nValues), label = "p = "+ ("%.2f"% p))
        savetxt("../results/1e/data-" + idString + ".dat", [sValues[:-1], nValues])
        print p
        
    xlabel(r"$\log_{10}[s]$")
    ylabel(r"$\log_{10}[n(s,p)]$")
    legend(prop={'size':11})
    xlim(xmin=1)
    grid()
    savefig("../results/1e/n-vs-s-" + idString + ".pdf")
    
    idString = "L" + str(Lx) +"x"+str(Ly) + "-nsamples" + str(nSamples)
    figure()



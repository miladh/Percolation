# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 15:59:17 2013

author: Milad H. Mobarhan 
"""

from pylab import *
import percolationLibrary as PL


pc = 0.59275
nSamples = 100
exponentList = range(4,12)
L = zeros(len(exponentList))
MValues  = zeros(len(exponentList))
PValues  = zeros(len(exponentList))
PiValues = zeros(len(exponentList))
idString = "nsamples" + str(nSamples)


kIndex=0
for k in exponentList:
    L[kIndex] = 2**k 
    kIndex+=1

Lx = L
Ly = L


for LIndex in range(len(L)):
    
    print Lx[LIndex]
    
    PiValues[LIndex], PValues[LIndex] = PL.clusterDensity(nSamples,Lx[LIndex],Ly[LIndex],pc)        
    MValues[LIndex] = PValues[LIndex]*Lx[LIndex]*Ly[LIndex]
    
    
plot(log10(L), log10(MValues))
#plot(pValues, PiValues, label  = "$\Pi(L,p)$, L = "+ ("%.2f"% Lx))
        
xlabel(r"$L$")
ylabel(r"$M(L)$")
grid()
savefig("../results/1h/M-vs-L-" + idString + ".pdf")

dMValues = diff(log10(MValues))/diff(log10(L))
figure()
plot(log10(L[:-1]), dMValues)      
xlabel(r"$L$")
ylabel(r"$\frac{dM(L)}{dL}$")
grid()
savefig("../results/1h/dM-vs-L-" + idString + ".pdf")

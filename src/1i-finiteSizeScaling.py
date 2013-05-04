# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 10:49:57 2013

@author: Milad H. Mobarhan
"""

from pylab import *
import percolationLibrary as PL
close("all")

npValues = 100
pValues  = linspace(0.54, 0.62, npValues)
nSamples = array([1000, 1000, 1000, 200, 200, 200]) 
#nSamples /= 100
L         = array([25,   50,   100, 200, 400, 800]) 

PValues  = zeros((len(L), npValues))
PiValues = zeros((len(L), npValues))
idString = "nsamples" + str(nSamples)


Lx = L
Ly = L

pPi0_8Values = zeros(len(L))
pPi0_3Values = zeros(len(L))

figure()
for LIndex in range(len(L)):
    print Lx[LIndex]
    
    for pIndex in range(len(pValues)):
        p = pValues[pIndex]    
        print p
        PiValues[LIndex,pIndex], PValues[LIndex,pIndex] = PL.clusterDensity(nSamples[LIndex],Lx[LIndex],Ly[LIndex],p)  

   
    
    lastIndex = where(PiValues[LIndex,:] <= 0.8)[0][-1]
    pPi0_8Values[LIndex] = pValues[lastIndex]
    
    lastIndex = where(PiValues[LIndex,:] <= 0.3)[0][-1]
    pPi0_3Values[LIndex]= pValues[lastIndex]
    
    plot(pValues, PiValues[LIndex])

figure()
plot(L, pPi0_3Values,'-o' ,label = "x = 0.3")
plot(L, pPi0_8Values,'-*', label = "x = 0.8")   
xlabel(r"$L$")
ylabel(r"$p_{\Pi = x}$")
grid()
legend()
savefig("../results/1i/pPi-vs-L-" + idString + ".pdf")

#1j) Findig my:
figure()
pPidiff = pPi0_8Values - pPi0_3Values
xlabel(r"$\log_{10}(L)$")
ylabel(r"$\log_{10}(p_{\Pi = 0.8}-p_{\Pi = 0.3})$")
grid()
plot(log10(L), log10(pPidiff),'-o')
savefig("../results/1i/dpPi-vs-L-" + idString + ".pdf")
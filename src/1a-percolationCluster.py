# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 17:42:41 2013

@author: Milad H. Mobarhan
"""
from pylab import *
import percolationLibrary as PL

close("all")


LList =array([100])

for L in LList:
    print "L = ", L
    Lx = L
    Ly = L
    
    npValues = 20
    nSamples = 1000
    pValues  = linspace(0.0, 1.0, npValues)
    PValues  = zeros(npValues)
    PiValues = zeros(npValues)
    
    for pIndex in range(len(pValues)):
        p = pValues[pIndex]
        PiValues[pIndex], PValues[pIndex] = PL.clusterDensity(nSamples,Lx,Ly,p)        
        print p
 
    idString = "L" + str(Lx) + "-nsamples" + str(nSamples)
    plot(pValues, PValues, label  = "$P(L,p)$, L = "+ ("%.2f"% Lx))
    plot(pValues, PiValues, label  = "$\Pi(L,p)$, L = "+ ("%.2f"% Lx))
    
xlabel(r"$p$")
ylabel(r"$P(L,p)$")
ylim(0,1.1)
xlim(pValues.min(),1.0)
grid()
legend(loc=2)
savefig("../results/1a/P-vs-p-" + idString + ".pdf")



#1b) Findig beta
pc = 0.59275
pNew    = zeros(npValues)
PNew    = zeros(npValues)

i=0
for p in pValues:
    if p > pc:
        pNew[i] = p
        i+=1
i=0       
for P in PValues:
    if P > 0:
      PNew[i] = P 
      i+=1

logpValues = log(pNew - pc)
logPValues = log(PNew)    
        
figure()
plot(logpValues, logPValues)
title("loglog for L = " + str(Lx))
xlabel(r"$\log(p - p_c)$")
ylabel(r"$\log(P(L,p))$")
grid()
savefig("../results/1a/P-vs-p-" + idString + "-loglog.pdf")
dlogPValues= diff(logPValues)/diff(logpValues)

figure()
plot(logpValues[1:], dlogPValues)
title("derivative of loglog for L = " + str(Lx))
xlabel(r"$\log(p - p_c)$")
ylabel(r"$\frac{d}{d(\log(p - p_c))} \log(P(L,p))$")
grid()

savefig("../results/1a/P-vs-p-" + idString + "-loglogderivative.pdf")

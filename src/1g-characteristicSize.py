# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 18:47:58 2013

author: Milad H. Mobarhan 
"""

from pylab import *
from scipy.ndimage import measurements
import percolationLibrary as PL

close("all")

Lx = 500
Ly = Lx
npValues = 5
nSamples = 1000
pc = 0.59275

Sxi = []
Pxi = []

idString = "L" + str(Lx) +"x"+ str(Ly) + "-nsamples" + str(nSamples)

pValues  = linspace(0.50, 0.58, npValues)

        
for pIndex in range(len(pValues)):
    p = pValues[pIndex]    
    
    sValues,nValues = PL.numberDensity(nSamples, Lx, Ly, p)
    sValuesPc,nValuesPc = PL.numberDensity(nSamples, Lx, Ly, pc)
    
    
    nRatio = nValues.astype(float)/nValuesPc
    cutoff = ones(len(nRatio))*-0.5

    for i in range(len(nRatio)):
        if cutoff[i]-0.1 <= log10(nRatio[i]) <= cutoff[i]+0.1:
            Sxi.append(sValues[i])
            Pxi.append(p) 
        
    
    plot(log10(sValues[:-1]),cutoff)    
    plot(log10(sValues[:-1]), log10(nRatio),'-*' ,label = "p = "+ ("%.2f"% p))
#    savetxt("../results/1e/data-" + idString + ".dat", [sValues[:-1], nValues])
    print p
   
   
xlabel(r"$\log_{10}[s]$")
ylabel(r"$\log_{10}[n(s,p)/n(s,p_c)]$")
legend(prop={'size':11})
xlim(xmin=1)
grid()
savefig("../results/1g/n-vs-s-" + idString + ".pdf")



SxiValues = array(Sxi)
PxiValues = array(Pxi)

figure()
plot(PxiValues,SxiValues,'-*')
ylabel(r"$S_{\xi}$")
xlabel(r"$p$")
grid()
savefig("../results/1g/sxi-vs-p-" + idString + ".pdf")



#Finding sigma:
logSxiValues = log10(SxiValues)
logPValues = log10(abs(PxiValues-pc))  
   

figure()
plot(logPValues, logSxiValues)
xlabel(r"$\log_{10}(|p - p_c|)$")
ylabel(r"$\log_{10}(s_\xi)$")
grid()
savefig("../results/1g/logsxi-vs-p-" + idString + ".pdf")   
   
   
dlogSxiValues = diff(logSxiValues)/diff(logPValues)
figure()
plot(logPValues[:-1], dlogSxiValues, '-*')
ylabel(r"$S_{\xi}$")
xlabel(r"$p$")
xlabel(r"$\log_{10}(|p - p_c|)$")
ylabel(r"$\frac{d}{d(\log_{10}(|p - p_c|))} \log_{10}(s_\xi)$")

grid()
savefig("../results/1g/dlogsxi-vs-p-" + idString + ".pdf")








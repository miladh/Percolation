# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 11:19:57 2013

@author: Svenn-Arne Dragly (small modification by Milad H. Mobarhan)
"""

from pylab import *
from scipy.ndimage import measurements
import percolationLibrary as PL
close("all")


makePlots  = False
pc = 0.59275
nSamples = 100
L = range(25,1000,200)
Lx = L
Ly = L

nSinglyConnctedSites = zeros(len(L))

for LIndex in range(len(L)):
    print "L X L = ", Lx[LIndex]

    for sample in range(nSamples):
        nCount = 0
        perc = []
        while (len(perc)==0):
            nCount = nCount + 1
            if (nCount >1000):
                print "Couldn't make percolation cluster..."
                break
            
            lattice = rand(Lx[LIndex],Ly[LIndex]) < pc
            labelMatrix, nClusters = measurements.label(lattice)
            perc_x = intersect1d(labelMatrix[0,:],labelMatrix[-1,:])
            perc = perc_x[where(perc_x > 0)]
            print nCount
        
        if len(perc) > 0:
            labelList = arange(nClusters + 1)
            areaList = measurements.sum(lattice, labelMatrix, index=labelList)
            areaLabelMatrix = areaList[labelMatrix]
            maxArea = areaList.max()
            spanningCluster = (labelMatrix == perc[0])
            
                   
            # Run walk on this cluster
            l,r = PL.walk(spanningCluster)
            
            singlyConnectedSites = l*r 
            nSinglyConnctedSites[LIndex] += size(where(singlyConnectedSites > 0)[0])
            zadd = spanningCluster + singlyConnectedSites
            
    nSinglyConnctedSites[LIndex]=nSinglyConnctedSites[LIndex]/nSamples
        
    if makePlots:
        figure()
        title("Spanning Cluster")
        imshow(spanningCluster, interpolation='nearest', origin='upper')
        savefig("../results/1l/spanningCluster.pdf") 

        figure()
        title("Left turning walker")        
        imshow(l, interpolation='nearest', origin='upper')
        
        figure()
        title("Right turning walker")
        imshow(r, interpolation='nearest', origin='upper')
        colorbar()
        
        figure()
        title("Singly connected sites")
        imshow(singlyConnectedSites, interpolation='nearest', origin ='upper')
        colorbar()
        
        figure()
        title("")
        imshow(zadd, interpolation='nearest', origin='upper')
        colorbar()
        show()
        

plot(log10(L), log10(nSinglyConnctedSites), '-o')
xlabel(r"$\log_{10}(L)$")
ylabel(r"$\log_{10}(M_{SC})$")
grid()
savefig("../results/1l/M-vs-L.pdf")
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 19:45:17 2013

author: Milad H. Mobarhan 
"""

from pylab import *
from scipy.ndimage import measurements
close("all")

Lx = 5
Ly = 5
r = rand(Lx,Ly)
p = 0.5
lattice = r < p

labelMatrix, nClusters = measurements.label(lattice)
labelList = arange(labelMatrix.max() + 1) 
shuffle(labelList) 
shuffledLabelMatrix = labelList[labelMatrix] 
imshow(shuffledLabelMatrix, origin='lower', interpolation='nearest') 
colorbar()
title("Percolation matrix")

# Calculate areas
figure()
areaList = measurements.sum(lattice, labelMatrix, index= arange(labelMatrix.max() + 1) )
areaLabelMatrix = areaList[labelMatrix]
areaImg = imshow(areaLabelMatrix, origin='lower', interpolation='nearest')
colorbar()
title("Clusters by area")
show()

# Bounding box
sliced = measurements.find_objects(areaLabelMatrix == areaLabelMatrix.max())
if(len(sliced) > 0):
    sliceX = sliced[0][1]    #cols are the x-values
    sliceY = sliced[0][0]    #rows are the y-values
    plotxlim=areaImg.axes.get_xlim()
    plotylim=areaImg.axes.get_ylim()
    plot([sliceX.start, sliceX.start, sliceX.stop, sliceX.stop, sliceX.start], \
         [sliceY.start, sliceY.stop, sliceY.stop, sliceY.start, sliceY.start], \
         color="red")
    xlim(plotxlim)
    ylim(plotylim)
show()

#idString = "L" + str(Lx) + "-nsamples" + str(nSamples)  
#savefig("../results/1a/P-vs-p-" + idString + ".pdf")

#figure()
#plot(logpVals, logPVals)
#title("loglog for L = " + str(Lx))
#xlabel(r"$\log(p - p_c)$")
#ylabel(r"$\log(P(L,p))$")
#grid()
#savefig("../results/1a/P-vs-p-" + idString + "-loglog.pdf")
#
#figure()
#dPdp = (logPVals[1:] - logPVals[:-1]) / (logpVals[1:] - logpVals[:-1])
#plot(logpVals[:-1], dPdp)
#title("derivative of loglog for L = " + str(Lx))
#xlabel(r"$\log(p - p_c)$")
#ylabel(r"$d\log(P)/d\log(p-p_c)$")
#grid()
#savefig("../results/1a/P-vs-p-" + idString + "-loglogderivative.pdf")
#
#savetxt("../results/1a/data-" + idString + ".dat", [pValues, PValues])
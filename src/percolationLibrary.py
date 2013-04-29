# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 20:26:41 2013

author: Milad H. Mobarhan 
"""
from pylab import *
from scipy.ndimage import measurements


def clusterDensity(nSamples, Lx, Ly, p):
    PValues  = 0
    PiValues = 0
    latticeArea = Lx*Ly
    for sample in range(nSamples):
            r = rand(Lx,Ly)
            lattice = r < p
            labelMatrix, nClusters = measurements.label(lattice)
            
            if(Lx==1):
                verticalPercolation = set()
            else:
                verticalPercolation   = set(labelMatrix[0,:]).intersection(labelMatrix[-1,:])
                
            if(Ly==1):
                horisontalPercolation = set()
            else:
                horisontalPercolation = set(labelMatrix[:,0]).intersection(labelMatrix[:,-1])
             
            if(Lx==1 and Ly==1):
                horisontalPercolation = set(labelMatrix[:,0]).intersection(labelMatrix[:,-1])
                verticalPercolation   = set(labelMatrix[0,:]).intersection(labelMatrix[-1,:])
                
                
            spanningClusterLabelList = array(list(horisontalPercolation|verticalPercolation))
            spanningClusterLabelList = delete(spanningClusterLabelList, where(spanningClusterLabelList == 0.0))
            
            if(len(spanningClusterLabelList) > 0.0):
                PiValues += 1.0
                areaList  = measurements.sum(lattice, labelMatrix, index = spanningClusterLabelList)               
                
                for area in areaList:                
                    PValues += area / latticeArea         
        
        
    PiValues /= nSamples
    PValues  /= nSamples
    
    return(PiValues, PValues) 



def numberDensity(nSamples, Lx, Ly, p):

    latticeArea = Lx*Ly
    Bins = logspace(0, log10(latticeArea))
    nBins = len(Bins)
    n = zeros(nBins - 1)    
    
    for sample in range(nSamples):
                r = rand(Lx,Ly)
                lattice = r < p
                labelMatrix, nClusters = measurements.label(lattice)
                labelList = arange(labelMatrix.max() + 1) 
                
                if(Lx==1):
                    verticalPercolation   = set()
                else:
                    verticalPercolation   = set(labelMatrix[0,:]) & set(labelMatrix[-1,:])
                    
                if(Ly==1):
                    horisontalPercolation = set()
                else:
                    horisontalPercolation = set(labelMatrix[:,0]) & set(labelMatrix[:,-1])
                 
                if(Lx==1 and Ly==1):
                    horisontalPercolation = set(labelMatrix[:,0]) & set(labelMatrix[:,-1])
                    verticalPercolation   = set(labelMatrix[0,:]) & set(labelMatrix[-1,:])
                    
                    
                spanningClusterLabelList = horisontalPercolation|verticalPercolation
                FiniteClusterLabelList = array(list(set(labelList) - spanningClusterLabelList))
                
                
                if(len(FiniteClusterLabelList) > 0.0):
                    areaList  = measurements.sum(lattice, labelMatrix, index = FiniteClusterLabelList)  
                else: 
                    areaList = zeros(len(FiniteClusterLabelList))
                
                
                
                data = histogram(areaList, bins = Bins)
                nValues = data[0]
                sValues = data[1]    
                
                ds = sValues[1:]-sValues[:-1]  
                nValues = nValues.astype(float)/(ds*len(areaList))           
                
                nValues = nValues.astype(float) /latticeArea
                n += nValues
                
    n = n/nSamples    
    
    return(Bins,n) 
    

def showMatrix(Lx,Ly,p): 
    r = rand(Lx,Ly)
    lattice = r < p
    
    labelMatrix, nClusters = measurements.label(lattice)
    labelList = arange(labelMatrix.max() + 1) 
    shuffle(labelList) 
    shuffledLabelMatrix = labelList[labelMatrix] 
    imshow(shuffledLabelMatrix, origin='lower', interpolation='nearest') 
    colorbar()
    title("Percolation matrix")
    
    return(lattice,labelMatrix) 
    

def clusterArea(Lx,Ly,p):

    lattice,labelMatrix= showMatrix(Lx,Ly,p)
    figure()
    areaList = measurements.sum(lattice, labelMatrix, index= arange(labelMatrix.max() + 1) )
    areaLabelMatrix = areaList[labelMatrix]
    imshow(areaLabelMatrix, origin='lower', interpolation='nearest')
    colorbar()
    title("Clusters by area")
    show()

    
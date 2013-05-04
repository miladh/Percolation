# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 11:19:57 2013

author: Svenn-Arne Dragly (small modification by Milad H. Mobarhan)
"""
from pylab import *
from scipy.ndimage import measurements
import percolationLibrary as PL
close("all")


makePlots  = False
pc = 0.59275
nSamples = 50
L = range(25,1000,200)
Lx = L
Ly = L
idstring = "L" + str(L) + "-p" + ("%.2f" % pc).replace(".", "_")

nSinglyConnctedSites = zeros(len(L))
nBackbone = zeros(len(L))
nDanglingEnds = zeros(len(L))

for LIndex in range(len(L)):
    print "L X L = ", Lx[LIndex]

    for sample in range(nSamples):
        nCount = 0
        perc = []   
        while (len(perc) == 0):
            nCount = nCount + 1
            if (nCount >100):
                print "Couldn't make percolation cluster..."
                break
            
            lattice = rand(Lx[LIndex],Ly[LIndex]) < pc
            labelMatrix, nClusters = measurements.label(lattice)
            perc_x = intersect1d(labelMatrix[0,:],labelMatrix[-1,:])
            perc = perc_x[where(perc_x > 0)]
            print "Percolation attempt", nCount
            
        if len(perc) > 0: 
            spanningCluster = asarray((labelMatrix == perc[0]))
            spanningClusterT = spanningCluster.T # Transpose
            
            # Generate bond lattice from this
            g = PL.siteToBond ( spanningClusterT )
            
            #Generate conductivity matrix
            p, c_eff = PL.FIND_COND (g, Lx[LIndex], Ly[LIndex])
            
            #Transform this onto a nx x ny lattice
            x = PL.colToMat ( p , Lx[LIndex] , Ly[LIndex] )
            P = x * spanningClusterT 
            
            g1 = g[:,0]
            g2 = g[: ,1]
            z1 = PL.colToMat( g1 , Lx[LIndex] , Ly[LIndex] )
            z2 = PL.colToMat( g2 , Lx[LIndex] , Ly[LIndex] )
            
            
            # Calculate flux from top to down 
            #(remember that flux is the negative of the pressure difference)
            f2 = zeros ( (Lx[LIndex] , Ly[LIndex] ))
            for iy in range(Ly[LIndex] -1):
                f2[: , iy ] = ( P [: , iy ] - P [: , iy +1]) * z2 [: , iy ]
            
            # Calculate flux from left to right
            #(remember that flux is the negative of the pressure difference)
            f1 = zeros ( (Lx[LIndex] , Ly[LIndex] ))
            for ix in range(Lx[LIndex]-1):
                f1[ ix ,:] = ( P [ ix ,:] - P [ ix +1 ,:]) * z1 [ ix ,:]
               
            #Find the sum of absolute fluxes in and out of each site
            fn = zeros (( Lx[LIndex] , Ly[LIndex] ))
            fn = fn + abs ( f1 )
            fn = fn + abs ( f2 )
            
            # Add for each column, except the leftmost one, the up-down flux, but offset
            fn [: ,1: Ly[LIndex] ] = fn [: ,1: Ly[LIndex] ] + abs ( f2 [: ,0: Ly[LIndex] -1])
            # For the left-most one, add the inverse pressure multiplied
            # with the spanning cluster bool information
            fn [: ,0] = fn [: ,0] + abs (( P [: ,0] - 1.0)*( spanningClusterT [: ,0]))
            # For each row except the topmost one, add the left-right flux, but offset
            fn [1: Lx[LIndex] ,:] = fn [1: Lx[LIndex] ,:] + abs ( f1 [0: Lx[LIndex] -1 ,:])
           
           
            singlyDiff = fn.max() - (sum(fn[:,0]) + sum(fn[:,-1]))
            if abs(singlyDiff) > 1e-8:
                singlyConnected = False
                scLimit = inf
            else:
                singlyConnected = True
                scLimit = fn.max()
            print "Has singly connected bonds:", singlyConnected
            
            zsc = (fn > scLimit - 1e-6)        
            zbb = (abs(fn) > 1e-6)
            zde = (spanningClusterT - zbb)
            
            
            nSinglyConnctedSites[LIndex] += size(where(zsc > 0)[0])
            nBackbone[LIndex] += size(where(zbb > 0)[0])
            nDanglingEnds[LIndex] += size(where(zde > 0)[0])
            
            
            #Plotting: 
            if makePlots:    
                figure()
                imshow(zsc * 3 + (zbb - zsc) * 2 + (spanningClusterT - zbb), interpolation='nearest')
                title ( "Singly connected, backbone and dead ends")
                savefig("../results/1o/sc-bb-de-" + idstring + ".pdf", dpi=300)
                grid(color="white")
                print sum(c_eff)        
                    
                    
                figure()
                ax = subplot(221)
                ax.set_adjustable('box-forced')
                imshow(spanningClusterT, interpolation='nearest')
                title("Spanning cluster")
                grid(color="white")
                
                ax3 = subplot(222, sharex=ax, sharey=ax)
                ax3.set_adjustable('box-forced')
                imshow(P, interpolation='nearest')
                title("Pressure")
                colorbar()
                grid(color="white")
                
                ax3 = subplot(223, sharex=ax, sharey=ax)
                ax3.set_adjustable('box-forced')
                imshow(fn, interpolation='nearest')
                title ( " Flux " )
                colorbar()
                grid(color="white")
                
                zfn = fn > fn.max() - 1e-6
                zbb = ( spanningClusterT + 2* zfn )
                zbb = zbb / zbb.max()
                ax3 = subplot(224, sharex=ax, sharey=ax)
                ax3.set_adjustable('box-forced')
                imshow(zbb, interpolation='nearest')
                title ( " BB and DE ")
                grid(color = "white")
                savefig("../results/1o/spanning-cluster-" + idstring + ".pdf", dpi=300)
                show()
          
          
          
nSinglyConnctedSites/=nSamples
nBackbone     /=nSamples
nDanglingEnds /=nSamples

plot(log10(L), log10(nSinglyConnctedSites),'-o', label = "Singly connected")
plot(log10(L), log10(nBackbone),'-o',label = "Backbone")
plot(log10(L), log10(nDanglingEnds),'-o', label = "Dangling ends")
xlabel(r"$\log_{10}(L)$")
ylabel(r"$\log_{10}(M_{SC})/\log_{10}(M_{DE})/\log_{10}(M_{BB})$")
grid()
legend(loc=2)
savefig("../results/1o/M-vs-L.pdf")
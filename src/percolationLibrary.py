# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 20:26:41 2013

author: Milad H. Mobarhan 
"""
from pylab import *
from scipy.ndimage import measurements
from scipy.sparse import spdiags, dia_matrix, coo_matrix
from scipy.sparse.linalg import spsolve

############################################################################
def FIND_COND (A , X , Y ):
# Sets up Kirchoff ’ s equations for the 2 D lattice A .
# A has X * Y rows and 2 columns . The rows indicate the site ,
# the first column the bond perpendicular to the flow direction
# and the second column the bond parallel to the flow direction .
# Calculates the effective flow conductance Ceff of the
# lattice A as well as the pressure P in every site .

    P_in = 1.
    P_out = 0.
    # Calls MK_EQSYSTEM .
    B,C = MK_EQSYSTEM (A , X , Y )

    # Kirchhoff ’ s equations solve for P
    P = spsolve(B, C)
    # The pressure at the external sites is added
    # ( Boundary conditions )
    P = concatenate((P_in * ones (X), P,  P_out * ones (X)))
    
    # Calculate Ceff
    Ceff = (P[-1-2*X+1:-1-X] - P_out).T * A[-1-2*X+1:-1-X, 1] / ( P_in - P_out )
    
    return P , Ceff
    
    
############################################################################
def MK_EQSYSTEM (A , X , Y ):
# The return values are [B , C ] where B * x = C . This is solved
# for the site pressure by x = B \ C .

    # Total no of internal lattice sites
    sites = X *( Y - 2)

    # Allocate space for the nonzero upper diagonals
    main_diag = zeros(sites)
    upper_diag1 = zeros(sites - 1)
    upper_diag2 = zeros(sites - X)
    
    # Calculates the nonzero upper diagonals
    main_diag = A[X:X*(Y-1), 0] + A[X:X*(Y-1), 1] + A[0:X*(Y-2), 1] + A[X-1:X*(Y-1)-1, 0]
    upper_diag1 = A [X:X*(Y-1)-1, 0]
    upper_diag2 = A [X:X*(Y-2), 1]
    main_diag[where(main_diag == 0)] = 1
    
    # Constructing B which is symmetric , lower = upper diagonals .
    B = dia_matrix ((sites , sites)) # B *u = t
    B = - spdiags ( upper_diag1 , -1 , sites , sites )
    B = B + - spdiags ( upper_diag2 ,-X , sites , sites )
    B = B + B.T + spdiags ( main_diag , 0 , sites , sites )
    
    # Constructing C
    C = zeros(sites)
    C[0:X] = A[0:X, 1]
    C[-1-X+1:-1] = 0*A [-1 -2*X + 1:-1-X, 1]
    
    return B , C
    
    
############################################################################
def siteToBond ( z ):
    # Function to convert the site network z (L , L ) 
    # into a ( L *L ,2) bond network
    # g [i,0] gives bond perpendicular to direction of flow
    # g [i,1] gives bond parallel to direction of flow
    # z [ nx , ny ] -> g [ nx * ny , 2]

    nx = size (z ,1 - 1)
    ny = size (z ,2 - 1)
    gg_r = zeros ((nx , ny))
    gg_d = zeros ((nx , ny ))
    
    gg_r [:, 0:ny - 1] = z [:, 0:ny - 1] * z [:, 1:ny]
    gg_r [: , ny  - 1] = z [: , ny  - 1]
    gg_d [0:nx - 1, :] = z [0:nx - 1, :] * z [1:nx, :]
    gg_d [nx - 1, :] = 0

    # Then , concatenate gg onto g
    g = zeros((nx *ny ,2))
    g [:, 0] = gg_d.reshape(-1,order='F').T
    g [:, 1] = gg_r.reshape(-1,order='F').T
    
    return g
    
############################################################################  
def colToMat (z, x, y):
    # Convert z ( x * y ) into a matrix of z (x , y )
    # Transform this onto a nx x ny lattice
    g = zeros ((x , y))

    for iy in range(1,y):
        i = (iy - 1) * x + 1
        ii = i + x - 1
        g[: , iy - 1] = z[ i - 1 : ii]
    
    return g

############################################################################
def walk(z):   
    # Returns: left & right: 
    # nr of times walker passes a site
    # First, ensure that array only has one contact point at left and
    # right : topmost points chosen
    nx = z.shape[0]
    ny = z.shape[1]
    i = where(z > 0)
    ix0 = 0 # starting row for walker is always 0
    iy0 = i[1][where(i[0] == 0)][0] # starting col (first element where there is a matching row which is zero)
    print "Starting walk in x=" + str(ix0) + " y=" + str(iy0)
    
    # Initilize directions
    directions = zeros((4,2), int)
    directions [0,0] = -1  # north
    directions [0,1] =  0
    directions [1,0] =  0  # west
    directions [1,1] = -1
    directions [2,0] =  1  # south
    directions [2,1] =  0
    directions [3,0] =  0  # east
    directions [3,1] =  1
    
    left = zeros((nx,ny),int)
    right = zeros((nx,ny),int)    
    
    # Left turning walker
    nwalk = 1
    ix = ix0
    iy = iy0
    direction = 0 # 0 = north, 1 = west, 2 = south, 3 = east
    while (nwalk > 0):
        left[ix,iy] = left[ix,iy] + 1
        # Turn left until you find an occupied site
        nfound = 0
        while (nfound == 0):
            direction = direction - 1
            if (direction < 0):
                direction = direction + 4
            
            # Check this direction
            iix = ix + directions[direction,0]
            iiy = iy + directions[direction,1]
            
            if (iix >= nx):
                nwalk  = 0 # Walker escaped
                nfound = 1
                iix = nx
                ix1 = ix
                iy1 = iy
                    
            # Is there a site here?
            elif(iix >= 0):
                if(iiy >= 0):
                    if (iiy < ny):
                        if (z[iix,iiy] > 0): # there is a site here, move here
                            ix = iix
                            iy = iiy
                            nfound = 1
                            direction = direction + 2    
                            if (direction > 3):
                                direction = direction - 4
        
    
    #Right turning walker
    nwalk = 1
    ix = ix0
    iy = iy0
    direction = 1 # 0 = north, 1 = west, 2 = south, 3 = east
    while(nwalk  > 0):
        right[ix,iy] = right[ix,iy] + 1
        # ix,iy
        # Turn right until you find an occupied site
        nfound = 0
        while (nfound==0):
            direction = direction + 1
            if (direction > 3):
                direction = direction - 4
            
            # Check this directionection
            iix = ix + directions[direction,0]
            iiy = iy + directions[direction,1]
            if (iix >= nx):
                if (iy >= iy1):
                    nwalk = 0 # Walker escaped
                    nfound = 1
                    iix = nx
                    
            # Is there a site here?
            elif(iix >= 0):
                if(iiy >= 0):
                    if (iiy < ny):
                        if (iix < nx):
                            if (z[iix,iiy]>0): # there is a site here, move here
                                ix = iix
                                iy = iiy
                                nfound = 1
                                direction = direction - 2
                                if (direction <0):
                                    direction = direction + 4
                                    
    return left, right



############################################################################
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


############################################################################
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
    
############################################################################
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
    
############################################################################
def clusterArea(Lx,Ly,p):
    lattice,labelMatrix= showMatrix(Lx,Ly,p)
    figure()
    areaList = measurements.sum(lattice, labelMatrix, index= arange(labelMatrix.max() + 1) )
    areaLabelMatrix = areaList[labelMatrix]
    imshow(areaLabelMatrix, origin='lower', interpolation='nearest')
    colorbar()
    title("Clusters by area")
    show()

    
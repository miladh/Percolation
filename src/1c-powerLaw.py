# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 10:37:42 2013

@author: Milad H. Mobarhan
"""

from pylab import*
close("all")
zValues = rand(1e6,1)**(-3+1);

Bins = logspace(log10(zValues.min()), log10(zValues.max()))
data = histogram(zValues, bins = Bins)

fz = data[0]
z  = data[1]

zdiff= z[1:]-z[:-1]

fzNormed =fz/ (zdiff*len(zValues))

figure()
plot(log10(z[:-1]),log10(fzNormed))
title("Plot of $f_Z(z)$")
xlabel(r"$z$")
ylabel(r"$f_Z(z)$")
grid()
savefig("../results/1c/distribution1.pdf")


Pz = cumsum(fz)
Pz.astype(float)
Pz=Pz/float(len(zValues))
figure()
plot(log10(z[:-1]),Pz)
title("Plot of $P (Z > z)$")
xlabel(r"$z$")
ylabel(r"$P (Z > z)$")
grid()  
savefig("../results/1c/cumulative-distribution.pdf")



dPz= diff(Pz)/diff(z[:-1])
figure()
plot(log10(z[:-2]),log10(dPz))
title("derivative of  $P (Z > z)$" )
xlabel(r"$\log_{10}(z)$")
ylabel(r"$\log_{10}(\frac{d}{dz} P (Z > z))$")
grid()
savefig("../results/1c/distribution2.pdf")



import scipy as sp
import numpy as np
import pylab as pl

D = np.loadtxt("gv-all.ascii")

p1 = D[::2]
p2 = D[1::2]

l = p1 - p2

Q = l[:,0:3] * l[:,3:6]

dl = l[1:] - l[:-1]

print p1
print p2

print l
print dl

dlx = dl[:,0:3]
dlv = dl[:,3:6]

L = dlx * dlv
print L
m =  np.mean(L[:,1])
print m
#L /= m
#L /= L[0,:]
#L -= 1

print L

pl.figure()
pl.plot(L[:,1])

pl.figure()
pl.plot(Q[:,1])

print L == Q

pl.show()


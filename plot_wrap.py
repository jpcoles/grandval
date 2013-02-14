import numpy as np
import pylab as pl

N = 10000
D = np.loadtxt("gv-all.ascii")

pl.ion()

xr = np.min(D[:,0]), np.max(D[:,0])
yr = np.min(D[:,3]), np.max(D[:,3])

for t in range(len(D) / N):
    X = D[t*N:(t+1)*N, 0]
    V = D[t*N:(t+1)*N, 3]
    pl.clf()
    pl.plot(X,V, ',')
    pl.xlim(xr)
    pl.ylim(yr)
    pl.draw()
    print "Writing frame", t
    pl.savefig("ps-%05i.png" % t)
    #raw_input()

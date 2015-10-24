from PageRankPipeline import *
import numpy as np
import matplotlib.pyplot as plt

SCALE=range(9,14)
EdgesPerVertex=16
Nfile=4
Niter=20

Nmax=pow(2*np.ones((len(SCALE))), SCALE)
M = EdgesPerVertex * Nmax

TimeArray=np.zeros((len(SCALE),4))

for i in xrange (0,len(SCALE)):
    
    TimeArray[i,0],TimeArray[i,1],TimeArray[i,2],TimeArray[i,3] = PageRankPipeline(SCALE[i], EdgesPerVertex, Nfile)

ElementMatrix = np.reshape(np.repeat(M, 4), (len(SCALE),4))

RateMatrix=ElementMatrix/TimeArray
Ratemat4 = RateMatrix[:,3]*Niter
RateMatrix[:,3] = Ratemat4

print RateMatrix
fig=plt.figure()
ax=fig.add_subplot(1,1,1)

line0, = ax.plot(M, RateMatrix[:,0], '--', linewidth=2)
line1, = ax.plot(M, RateMatrix[:,1], '--', linewidth=2)
line2, = ax.plot(M, RateMatrix[:,2], '--', linewidth=2)
line3, = ax.plot(M, RateMatrix[:,3], '--', linewidth=2)

plt.legend(['K0 rate', 'K1 rate', 'K2 rate', 'K3 rate'])
plt.xlabel('Number of Edges')
plt.ylabel('edges/second')

ax.set_xscale('log')
ax.set_yscale('log')

plt.show()

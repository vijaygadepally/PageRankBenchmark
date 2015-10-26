PLOT=0
from PageRankPipeline import *
import numpy as np

if PLOT:
    import matplotlib.pyplot as plt

SCALE=range(10,23)
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

f=open('python.out','w')
strtowrite="SCALE \t K0-edges-per-s \t K1-edges-per-s \t K2-edges-per-s \t K3-edges-per-s\n"
f.write(strtowrite)

#Write Ratematrix
out=""
for row in xrange (0,len(SCALE)):
    tmp="%d" % SCALE[row] + "\t" + "%0f" % RateMatrix[row,0] + "\t" + "%0f" % RateMatrix[row,1] + "\t" + "%0f" % RateMatrix[row,2]+ "\t" + "%0f" % RateMatrix[row,3]+chr(10)
    out+=tmp

f.write(out)

print RateMatrix

if PLOT:
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
    plt.savefig('python_pgresults.pdf')
    #plt.show()

########################################################
# PageRank Pipeline Benchmark
# Developer: Dr. Vijay Gadepally (vijayg@mit.edu)
# MIT
########################################################
# (c) <2015> Vijay Gadepally
########################################################

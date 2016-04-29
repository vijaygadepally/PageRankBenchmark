import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import csc_matrix
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigs
import time
from numpy import linalg as la
#Use the pandas package if available
#import pandas as pd

import scipy.io as sio


###################################################
###################################################
def KronGraph500NoPerm (SCALE, EdgesPerVertex):
    N=pow(2,SCALE)
    M = EdgesPerVertex * N

    A=0.57
    B=0.19
    C=0.19
    D=1-(A+B+C)

    ij=np.zeros((2,M))

    ab=A+B
    c_norm=C/(1-(A+B))
    a_norm=A/(A+B)

    for ib in range (0,SCALE):
        ii_num=np.random.uniform(0,1,M)
        ii_bit=ii_num>ab
        jj_num=np.random.uniform(0,1,M)
        jj_bit = jj_num > (c_norm * ii_bit + a_norm*np.invert(ii_bit))
        ijcomb = (pow(2,ib))*np.vstack((ii_bit,jj_bit)) #allow vertices from 0
        ij = ij+ijcomb

    StartVertex, EndVertex = np.vsplit(ij, 2)
    return StartVertex, EndVertex;

###################################################
###################################################
def StrArrayWrite(nparray, filename):

    f=open(filename,"wt",buffering=20*(1024**2))
    f=open(filename,"w")
    data = [str(float(row[0])) + '\t' + str(float(row[1])) + '\n' for row in nparray]
    f.write(''.join(data))
    f.close()

    #Use Pandas if you have it
    #pd.DataFrame(nparray).to_csv(filename, sep='\t', header=False, index=False)



###################################################
###################################################
def StrArrayRead(filename):

    f=open(filename,'r')
    edgelist = []
    with open(filename, 'r') as f:
        for line in f:
            edgelist.append(list(map(float, line.split('\t'))))
    f.close()
    return np.asarray(edgelist)

    #Use Pandas if you have it
    #return pd.read_csv(filename, delimiter='\t').as_matrix()


###################################################
###################################################
#@profile
def PageRankPipeline (SCALE, EdgesPerVertex, Nfile):

    Nmax=pow(2,SCALE)
    M = EdgesPerVertex * Nmax
    c=0.85
    Niter=20

    print "Number of Edges " + str(M) + ", Max Possible Vertex: " + str(Nmax)

    ###################################################
    ###################################################

    #Kernel 0: Generate random graph and save to file
    print "Kernel 0: Generate Graph, Write Edges"

    startTime=time.clock()

    for i in xrange (0,Nfile):
        fname= "data/K0/" + str(i) + ".tsv"
        print "   Writing: " + fname
        np.random.seed(i)
        u, v = KronGraph500NoPerm(SCALE,EdgesPerVertex/Nfile)
        uv=np.vstack((u,v))
        StrArrayWrite(np.transpose(uv), fname)

    K0time = time.clock() - startTime
    print "K0time " + str(K0time) + ", Edges/sec: " + str( M/K0time )

    ###################################################
    ###################################################

    #Kernel 1: Read, Sort, Write Edges
    print "Kernel 1: Read, Sort, Write Edges"
    startTime=time.clock()

    edgelist=np.empty((0,2))

    #Read into a single array
    for i in xrange (0,Nfile):
        fname= "data/K0/" + str(i) + ".tsv"
        print "   Reading:" + fname
        tmp = StrArrayRead(fname)
        edgelist=np.concatenate((edgelist,tmp),axis=0)

    #Sort by start edge
    u,v = np.hsplit(edgelist,2)
    ind=np.argsort(np.transpose(u))
    u= np.take(u, ind)
    v= np.take(v, ind)

    #Write Data to files
    j=1

    for i in xrange (0,Nfile):
        jEdgeStart = ((j-1) * u.size / Nfile)
        jEdgeEnd = ((j) * u.size/ Nfile)
        ind = range(jEdgeStart, jEdgeEnd)
        uu = np.take(u, ind)
        vv = np.take(v, ind)
        uuvv=np.vstack((uu,vv))
        fname= "data/K1/" + str(i) + ".tsv"
        print "   Writing: " + fname
        StrArrayWrite(np.transpose(uuvv), fname)
        j=j+1

    K1time = time.clock()-startTime
    print "K1time " + str(K1time) + ", Edges/sec: " + str( M/K1time )

    ###################################################
    ###################################################

    #Kernel 2: Read, Filter Edges
    print "Kernel 2: Read, Filter Edges"
    startTime=time.clock()

    edgelist=np.empty((0,2))

    #Read into a single array
    for i in xrange (0,Nfile):
        fname= "data/K1/" + str(i) + ".tsv"
        print "   Reading:" + fname
        tmp = StrArrayRead(fname)
        edgelist=np.concatenate((edgelist,tmp),axis=0)

    #Construct adjacency matrix
    u,v = np.hsplit(edgelist,2)
    d=np.ones(u.size).reshape(u.size,1)
    ut=np.squeeze(u)
    vt=np.squeeze(v)
    dt=np.squeeze(d)

    A=csc_matrix((dt, (ut,vt)), shape=(Nmax, Nmax))

    #Filter and weight data
    din = A.sum(axis=0)

    for row in np.ravel(np.nonzero(din==np.max(din))*1):
        A.data[A.indptr[row]:A.indptr[row+1]] =0

    for row in np.ravel(np.nonzero(din==1)*1):
        A.data[A.indptr[row]:A.indptr[row+1]] =0

    dout = A.sum(axis=1) #dout =sum(np.transpose(A))
    dinv=np.squeeze(np.asarray(1/dout))
    dinv[np.isinf(dinv)]=0
    dind=np.asarray(np.where(dinv>0))
    dval=dinv[dinv>0]

    D=csc_matrix((dval, (np.squeeze(dind),np.squeeze(dind))), shape=(Nmax,Nmax))

    #Perform filtering
    A= D*A

    K2time=time.clock()-startTime
    print "K2time " + str(K2time) + ", Edges/sec: " + str( M/K2time )

    #raise NameError('Die!')
    ###################################################
    ###################################################

    #Kernel 3: Compute PageRank.
    print "Kernel 3: Compute PageRank."
    startTime=time.clock()

    r=np.transpose(np.random.uniform(0,1,(Nmax,1)))
    r=r/la.norm(r,1)
    a= (1-c)/Nmax

    #r = ((c .* r) *A) + ((1-c) .* sum(r,2))
    for i in xrange (0,Niter):
        r = ( (c*r) * A ) + (a * r.sum(axis=1))
    
    r=r/la.norm(np.transpose(r),1)

    K3time=time.clock()-startTime
    print "Pagerank Sum= " + str(r.sum(axis=1))
    print "K3time " + str(K3time) + ", Edges/sec: " + str( Niter*M/K3time )

    #Check correctness
    #Check with eigenvector.
    #MATLAB output: sio.savemat('A.mat', {'vecs':vecs})
    #vals, vecs=eigs(c*np.transpose(A.todense()) + (1-c)/Nmax,k=1)
    #vecs=vecs/la.norm(vecs,1)
    #print np.max(abs(np.transpose(vecs)-r))
    
    return K0time, K1time, K2time, K3time;


########################################################
# PageRank Pipeline Benchmark
# Developer: Dr. Vijay Gadepally (vijayg@mit.edu)
# MIT
########################################################
# (c) <2015> Vijay Gadepally
########################################################

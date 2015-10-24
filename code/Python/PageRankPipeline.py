import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix
import time
from numpy import linalg as la


# import code
# import sys
# def keyboard(banner=None):
#     ''' Function that mimics the matlab keyboard command '''
#     # use exception trick to pick up the current frame
#     try:
#         raise None
#     except:
#         frame = sys.exc_info()[2].tb_frame.f_back
#     print "# Use quit() to exit :) Happy debugging!"
#     # evaluate commands in current namespace
#     namespace = frame.f_globals.copy()
#     namespace.update(frame.f_locals)
#     try:
#         code.interact(banner=banner, local=namespace)
#     except SystemExit:
#         return

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

    np.savetxt(filename, nparray, fmt='%16d', delimiter='\t', newline='\n')

    # fo=open(filename, "w")
    # nparray.tofile(fo, sep="\t", format="%s")
    # fo.close

###################################################
###################################################
def StrArrayRead(filename):

    #edgelist=np.fromfile(filename, sep='\t', dtype=np.float)
    #return edgelist

    # FASTER READ
    edgelist = []
    with open(filename, 'r') as f:
        for line in f:
            edgelist.append(list(map(float, line.split('\t'))))
    return np.asarray(edgelist)

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
    #edgelist=np.array([])

    #Read into a single array
    for i in xrange (0,Nfile):
        fname= "data/K0/" + str(i) + ".tsv"
        print "   Reading:" + fname
        tmp = StrArrayRead(fname)
        edgelist=np.concatenate((edgelist,tmp),axis=0)
        #edgelist = np.hstack((edgelist, tmp))

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

    #edgelist=np.array([])
    edgelist=np.empty((0,2))

    #Read into a single array
    for i in xrange (0,Nfile):
        fname= "data/K1/" + str(i) + ".tsv"
        print "   Reading:" + fname
        tmp = StrArrayRead(fname)
        edgelist=np.concatenate((edgelist,tmp),axis=0)
        #edgelist = np.hstack((edgelist, tmp))

    #Construct adjacency matrix
    u,v = np.hsplit(edgelist,2)
    d=np.ones(u.size).reshape(u.size,1)
    ut=np.squeeze(u)
    vt=np.squeeze(v)
    dt=np.squeeze(d)

    A=csr_matrix((dt, (ut,vt)), shape=(Nmax, Nmax))

    #Filter and weight data
    din = A.sum(axis=0)
    A[:,np.ravel(din==np.max(din))]=0 #A[:,din==max(din)]=0
    A[:,np.ravel(din==1)]=0 #A[:,din==1]=0
    dout = A.sum(axis=1) #dout =sum(np.transpose(A))
    dinv=np.squeeze(np.asarray(1/dout))
    dinv[np.isinf(dinv)]=0
    dind=np.asarray(np.where(dinv>0))
    dval=dinv[dinv>0]

    D=csr_matrix((dval, (np.squeeze(dind),np.squeeze(dind))), shape=(Nmax,Nmax))
    #Dinv=np.diag(dinv, 0)

    A = A.dot(D)

    K2time=time.clock()-startTime
    print "K2time " + str(K2time) + ", Edges/sec: " + str( M/K2time )

    #raise NameError('Die!')

    ###################################################
    ###################################################

    #Kernel 3: Compute PageRank.
    print "Kernel 3: Compute PageRank."
    startTime=time.clock()

    r=np.random.uniform(0,1,(Nmax,1))
    r=r/la.norm(r)

    a=(np.ones((Nmax,1)) * (1-c))/Nmax

    for i in xrange (0,Niter):
        r = A.dot(r)*c + a

    K3time=time.clock()-startTime
    print "Pagerank Sum= " + str(r.sum(axis=0))
    print "K3time " + str(K3time) + ", Edges/sec: " + str( M/K3time )

    return K0time, K1time, K2time, K3time;

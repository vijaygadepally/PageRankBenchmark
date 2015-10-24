import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix
import time
from numpy import linalg as la

###################################################
###################################################
#@profile
def KronGraph500NoPerm (SCALE, EdgesPerVertex):	
	N=pow(2,SCALE)
	M = EdgesPerVertex * N

	A=0.57
	B=0.19
	C=0.19
	D=1-(A+B+C)

	ij=np.ones((2,M))
	ab=A+B
	c_norm=C/(1-(A+B))
	a_norm=A/(A+B)

	for ib in xrange (0,SCALE):
		ii_num=np.random.uniform(0,1,M)
		ii_bit=ii_num>ab
		jj_num=np.random.uniform(0,1,M)
		jj_bit = jj_num > (c_norm * ii_bit + a_norm*np.invert(ii_bit))
		ijcomb = pow(2,ib)*np.vstack((ii_bit,jj_bit)) #allow vertices from 0
        ij = ij+ijcomb
        print ij
        StartVertex, EndVertex = np.vsplit(ij, 2)
	print StartVertex
	return StartVertex, EndVertex;

###################################################
###################################################
#@profile
def StrArrayWrite(nparray, filename):
	fo=open(filename, "w")
	nparray.tofile(fo, sep="\t", format="%s")
	#print nparray
	fo.close

###################################################
###################################################
#@profile
def StrArrayRead(filename):
	#fo=open(filename, "r")
	edgelist=np.fromfile(filename, sep='\t', dtype=np.float)
    #edgelist=np.load(filename)
	#edgelist = genfromtxt(filename, delimiter='')
	#edgelist=pd.read_csv(filename,delimiter='').values
	#fo.close
	return edgelist

###################################################
###################################################
#@profile
def PageRankPipeline (SCALE, EdgesPerVertex, Nfile):

#    SCALE=12
#    EdgesPerVertex=12
#    Nfile=4

    Nmax=pow(2,SCALE)
    M = EdgesPerVertex * Nmax
    c=0.15
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
        np.random.seed
        u, v = KronGraph500NoPerm(SCALE,EdgesPerVertex/Nfile)
        uv=np.vstack((u,v))
        StrArrayWrite(np.transpose(uv) , fname)

    K0time = time.clock() - startTime
    print "K0time " + str(K0time) + ", Edges/sec: " + str( M/K0time )


    ###################################################
    ###################################################

    #Kernel 1: Read, Sort, Write Edges
    print "Kernel 1: Read, Sort, Write Edges"
    startTime=time.clock()

    edgelist=np.array([])

    #Read into a single array
    for i in xrange (0,Nfile):
        fname= "data/K0/" + str(i) + ".tsv"
        print "   Reading:" + fname
        tmp = StrArrayRead(fname)
        edgelist=np.concatenate((edgelist,tmp),axis=0)
        #edgelist = np.hstack((edgelist, tmp))

    edgelist=edgelist.reshape((Nmax*EdgesPerVertex,2))

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

    #raise NameError('Die!')
    ###################################################
    ###################################################

    #Kernel 2: Read, Filter Edges
    print "Kernel 2: Read, Filter Edges"
    startTime=time.clock()

    edgelist=np.array([])       

    #Read into a single array
    for i in xrange (0,Nfile):
        fname= "data/K1/" + str(i) + ".tsv"
        print "   Reading:" + fname
        tmp = StrArrayRead(fname)
        edgelist=np.concatenate((edgelist,tmp),axis=0)
        #edgelist = np.hstack((edgelist, tmp))

    edgelist=edgelist.reshape((Nmax*EdgesPerVertex,2))

    #Construct adjacency matrix
    u,v = np.hsplit(edgelist,2)
    d=np.ones(u.size).reshape(u.size,1)
    ut=np.squeeze(u)
    vt=np.squeeze(v)
    dt=np.squeeze(d)

    print max(ut)
    print min(ut)
    A=csr_matrix((dt, (ut,vt)), shape=(Nmax, Nmax)).toarray()


    print "sparsity: " + str(np.count_nonzero(A)) + ' ' + str(Nmax)
    #Filter and weight data
    din = A.sum(axis=0)
    A[:,din==max(din)]=0
    A[:,din==1]=0
    dout = A.sum(axis=1)
    #dout =sum(np.transpose(A))
    dinv=1/dout
    dinv[np.isinf(dinv)]=0
    Dinv=np.diag(dinv, 0)
    A = A*Dinv
    #Dinv= lil_matrix((Nmax,Nmax))
    #Dinv.setdiag(dinv)
    #Dinv=np.reshape(np.repeat(dinv, Nmax), (Nmax,Nmax))
    #A = np.dot(Dinv, A)

    K2time=time.clock()-startTime
    print "K2time " + str(K2time) + ", Edges/sec: " + str( M/K2time )
    raise NameError('Die!')


    ###################################################
    ###################################################

    #Kernel 3: Compute PageRank.
    print "Kernel 3: Compute PageRank."
    startTime=time.clock()

    r=np.random.uniform(0,1,Nmax)
    r=r/la.norm(r)

    a=(np.ones(Nmax) * (1-c))/Nmax

    for i in xrange (0,Niter):
        r= A*(r*c) + a

    K3time=time.clock()-startTime
    print "K3time " + str(K3time) + ", Edges/sec: " + str( M/K3time )

    return K0time, K1time, K2time, K3time;

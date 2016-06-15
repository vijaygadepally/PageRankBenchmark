#
include("KronGraph500NoPerm.jl")
include("StrFileWrite.jl")
include("StrFileRead.jl")

function PageRankPipeline(SCALE,EdgesPerVertex,Nfile);

  Nmax = 2.^SCALE;                           # Max vertex ID.
  M = EdgesPerVertex .* Nmax;                # Total number of edges.
  myFiles = collect(1:Nfile).';              # Set list of files.
  #
  # Julia parallel version
  # Figure it out later: how to distribute the load
  # myFiles = global_ind(zeros(Nfile,1,map([Np 1],{},0:Np-1)));   # PARALLEL.
  tab = Char(9)
  nl = Char(10)
  Niter = 20                                      # Number of PageRank iterations.
  c = 0.15                                        # PageRank damping factor.

  println("Number of Edges: " * string(M) * ", Maximum Possible Vertex: " * string(Nmax));


  ########################################################
  # Kernel 0: Generate a Graph500 Kronecker graph and save to data files.
  ########################################################
  println("Kernel 0: Generate Graph, Write Edges");
  tic();
    for i in myFiles
      fname = "data/K0/" * string(i) * ".tsv";
      println("  Writing: " * fname);                          # Read filename.
      srand(i);                                                # Set random seed to be unique for this file.
      u, v = KronGraph500NoPerm(SCALE,EdgesPerVertex./Nfile);  # Generate data.

      writeuv(fname, u, v)
    end
  K0time = toq();
  println("K0 Time: " * string(K0time) * ", Edges/sec: " * string(M./K0time));


  ########################################################
  # Kernel 1: Read data, sort data, and save to files.
  ########################################################
  println("Kernel 1: Read, Sort, Write Edges");
  tic();

    # Read in all the files into one array.
    for i in myFiles
      fname = "data/K0/" * string(i) * ".tsv";
      println("  Reading: " * fname);  # Read filename.
      ut,vt = StrFileRead(fname);
      # Concatenate to u,v
      if i == 1
         u = ut; v = vt;
      else
         append!(u, ut)
         append!(v, vt)
      end
    end

    sortIndex = sortperm(u)                      # Sort starting vertices.
    u = u[sortIndex]                                  # Get starting vertices.
    v = v[sortIndex]                                  # Get ending vertices.

  K1time1 = toq();
  tic();
    # Write all the data to files.
    j = 1;                                                         # Initialize file counter.
    c = size(u,1)/length(myFiles)        # Compute first edge of file.
    for i in myFiles
      jEdgeStart = round(Int, (j-1)*c+1)# Compute first edge of file.
      jEdgeEnd = round(Int, j*c)          # Compute last edge of file.
      uu = sub(u,jEdgeStart:jEdgeEnd)                                 # Select start vertices.
      vv = sub(v,jEdgeStart:jEdgeEnd)                                 # Select end vertices.
      fname = "data/K1/" * string(i) * ".tsv"
      println("  Writing: " * fname)                              # Create filename.

      writeuv(fname, uu, vv)

      j = j + 1                                                   # Increment file counter.
    end

  K1time2 = toq();
  K1time = K1time1 + K1time2;
  println("K1 Time (reading):" * string(K1time1) * ", Edges/sec: " * string(M./K1time1));
  println("K1 Time (writing):" * string(K1time2) * ", Edges/sec: " * string(M./K1time1));
  println("K1 Time: " * string(K1time) * ", Edges/sec: " * string(M./K1time));


  ########################################################
  # Kernel 2: Read data, filter data.
  ########################################################
  println("Kernel 2: Read, Filter Edges");
  tic();
    # Read in all the files into one array.
    for i in myFiles
      fname = "data/K1/" * string(i) * ".tsv";
      println("  Reading: " * fname);                # Read filename.
      ut,vt = StrFileRead(fname);
      if i == 1
         u = ut; v = vt;                             # Initialize starting and ending vertices.
      else
         append!(u, ut)
         append!(v, vt)
         # Get the rest of starting and ending vertices.
      end
    end

    # Construct adjacency matrix.
    A = sparse(vec(u),vec(v),1.0,Nmax,Nmax)      # Create adjacency matrix.

    # Filter and weight the adjacency matrix.
    din = sum(A,1)                               # Compute in degree.
    A[find(din == maximum(din))]=0               # Eliminate the super-node.
    A[find(din == 1)]=0                          # Eliminate the leaf-node.
    dout = sum(A,2)                              # Compute the out degree.
    i = find(dout)                               # Find vertices with outgoing edges (dout > 0).
    DoutInvD = zeros(Nmax)        # Create diagonal weight matrix.
    DoutInvD[i] = 1./dout[i]
    scale!(DoutInvD, A)           # Apply weight matrix.
  K2time = toq();
  println("K2 Time: " * string(K2time) * ", Edges/sec: " * string(M./K2time));


  ########################################################
  # Kernel 3: Compute PageRank.
  ########################################################
  println("Kernel 3: PageRank");
  tic();

    r = rand(1,Nmax);                     # Generate a random starting rank.
    r = r ./ norm(r,1);                   # Normalize
    a = (1-c) ./ Nmax;                    # Create damping vector

    for i=1:Niter
        s = r * A
        scale!(s, c)
        r = s .+ (a * sum(r,2));                # Compute PageRank.
    end

    r=r./norm(r,1);

  K3time = toq();
  println("  Sum of PageRank: " * string(sum(r)) );     # Force all computations to occur.
  println("K3 Time: " * string(K3time) * ", Edges/sec: " * string(Niter.*M./K3time));

  return K0time,K1time,K2time,K3time

end

########################################################
# PageRank Pipeline Benchmark
# Architect: Dr. Jeremy Kepner (kepner@ll.mit.edu)
# Julia Translation: Dr. Chansup Byun (cbyun@ll.mit.edu)
# MIT
########################################################
# (c) <2015> Massachusetts Institute of Technology
########################################################

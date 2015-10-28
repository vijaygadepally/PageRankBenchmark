#!/usr/local/Julia/latest/bin/julia
#
using PyPlot

include("PageRankPipeline.jl")

# Workaround for Matlab bsxfun
function bsxfun_rdivide(M,Ktime)
   broadcast(./,M,Ktime)
end

SCALE = [10:16]';           # Scale of problem (transpose to match Matlab style).
EdgesPerVertex = 16;        # Average degree of each vertex (power of 2).
Nfile = 4;                  # Number of files to use (any power of 2).
Niter = 20;                 # Number of PageRank iterations.

Nmax = 2.^SCALE;            # Max vertex ID.
M = EdgesPerVertex .* Nmax; # Total number of edges.

# Run the PageRank Pipeline benchmark at each scale.
Ktime = zeros(4,length(SCALE));
for i=1:length(SCALE)
  Ktime[1,i],Ktime[2,i],Ktime[3,i],Ktime[4,i] = PageRankPipeline(SCALE[i],EdgesPerVertex,Nfile);
end

# Compute the rate.
Krate = bsxfun_rdivide(M,Ktime);
Krate[4,:] = Niter .* Krate[4,:];

# Figure out the plot later
# figure;
# loglog(M,Krate);  xlabel('number of edges');  ylabel('edges/second');
# legend('K0 Generate','K1 Sort','K2 Filter','K3 PageRank','Location','NorthEastOutside');


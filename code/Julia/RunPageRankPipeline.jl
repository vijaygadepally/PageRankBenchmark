#!/usr/local/Julia/latest/bin/julia
#
PLOTRESULT=0
if PLOTRESULT != 0
   using Plots
end

include("PageRankPipeline.jl")

# Workaround for Matlab bsxfun
function bsxfun_rdivide(M,Ktime)
   broadcast(./,M,Ktime)
end

SCALE = collect(10:16)';           # Scale of problem (transpose to match Matlab style).
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

if PLOTRESULT == 0
   # Write the result to a file
   # scale k0-edges-per-sec k1-edges-per-sec k2-edges-per-sec k3-edges-per-sec, tab separated
   header = ["scale", "k0-edges-per-sec", "k1-edges-per-sec", "k2-edges-per-sec", "k3-edges-per-sec" ]
   fid = open("julia.dat","w");
   writedlm(fid,header.','\t'); 
   writedlm(fid,[M; Krate].','\t'); 
   close(fid);
else
   # Julia plot uisng gadfly() aloing with Plots package
   # Requires web browser to display the plot
   # Works with Julia 0.4.0, Need to install the following packages and their dependencies
   #
   #  "Gadfly"            => v"0.3.18"
   #  "Plots"             => v"0.4.0"
   #
   gadfly()
   # styles = setdiff(supportedStyles(),[:auto])';
   labels = ["K0 Generate", "K1 Sort", "K2 Filter", "K3 PageRank"]';
   plot(M',Krate',style=:auto,label=map(string,labels),w=5)
   xaxis!("number of edges",:log10)
   yaxis!("edges/second",:log10)
end

########################################################
# PageRank Pipeline Benchmark
# Architect: Dr. Jeremy Kepner (kepner@ll.mit.edu)
# Julia Translation: Dr. Chansup Byun (cbyun@ll.mit.edu)
# MIT
########################################################
# (c) <2015> Massachusetts Institute of Technology
########################################################

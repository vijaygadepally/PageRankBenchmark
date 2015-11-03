DRAFT - PageRank Pipeline Benchmark  - DRAFT

1. Contents

PageRankPipeline/README.txt (this file)

PageRankPipeline/doc (documentation)
PageRankPipeline/doc/151012-PageRankPipelineBenchmark.pdf (description of benchmark)
PageRankPipeline/doc/SerialPerformanceMatlab.pdf (plot of serial matlab performance of benchmark)

PageRankPipeline/code (source code)
PageRankPipeline/code/Matlab (Matlab source code)
PageRankPipeline/code/Octave (Octave source code)
PageRankPipeline/code/Julia (Julia source code)
PageRankPipeline/code/Python (Python source code)
PageRankPipeline/code/R (R source code)
PageRankPipeline/code/Java (Java source code)


2. Execution

2.1 Matlab

Go to PageRankPipeline/code/Matlab directory.

Edit PageRankPipeline/code/Matlab/RunPageRankPipeline.m set SCALE and Nfile variables as desired.

Start matlab

type RunPageRankPipeline


2.2 Octave

Go to PageRankPipeline/code/Octave directory.

Edit PageRankPipeline/code/Octave/RunPageRankPipeline.m set SCALE and Nfile variables as desired.

Start Octave

type RunPageRankPipeline

2.3 Julia
 
Prerequisite: Julia 0.4.0 release
 
Go to PageRankPipeline/code/Julia directory.
 
Edit PageRankPipeline/code/Octave/RunPageRankPipeline.jl set SCALE and Nfile variables as desired.
 
Make sure Julia 0.4.0 version is used.
 
julia RunPageRankPipeline.jl


2.4 Python

python runPageRankPipeline.py


  You will need to install:  numpy scipy python-matplotlib
  E.g., in Fedora 22 do
     sudo dnf install numpy scipy python-matplotlib
 
2.4 R


2.5 Java

2.6 C++

 $ cd C++
 $ make
        (or "make CXX=icc" if you have the Intel compiler installed)
 $ ./runpagerankpipeline
 

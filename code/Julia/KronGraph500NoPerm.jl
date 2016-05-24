function KronGraph500NoPerm(SCALE,EdgesPerVertex)
#
#Graph500NoPerm: Generates graph edges using the same 2x2 Kronecker algorithm (R-MAT) as the Graph500 benchmark,
#                but no permutation of vertex labels is performed.
#IO user function.
#  Usage:
#    [StartVertex EndVertex] = Graph500NoPerm(SCALE,edgefactor)
#  Inputs:
#    SCALE = integer scale factor that sets the max number of vertices to 2^SCALE
#    EdgesPerVertex = sets the total number of edges to M = K*N;
# Outputs:
#    StartVertex = Mx1 vector of integer start vertices in the range [1,N]
#    EndVertex = Mx1 vector of integer end vertices in the range [1,N]

  N = 2.^SCALE                       # Set  power of number of vertices..

  M = round(Int, EdgesPerVertex .* N)     # Compute total number of edges to generate.

  A = 0.57; B = 0.19;  C = 0.19;   D = 1-(A+B+C)  # Set R-MAT (2x2 Kronecker) coefficients.

  ij = ones(Int, 2, M)         # Initialize index arrays.
  ab = A + B                 # Normalize coefficients.
  c_norm = C/(1 - (A + B))
  a_norm = A/(A + B)

  for j = 1:M
      for ib = 1:SCALE            # Loop over each scale.
          k = 1 << (ib-1)
          if rand() > ab
              ij[1,j] += k
              if rand() > c_norm
                  ij[2,j] += k
              end
          elseif rand() > a_norm
              ij[2,j] += k
          end
      end
  end

  StartVertex = sub(ij,1,:)     # Copy to output. (row vector)
  EndVertex =   sub(ij,2,:)       # Copy to output. (row vector)

  return StartVertex,EndVertex
end

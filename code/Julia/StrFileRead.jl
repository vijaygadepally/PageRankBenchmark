#
# Read Edges from file
#
function StrFileRead(fname)
    uv = readdlm(fname, ' ', Int)::Matrix{Int}
    n = size(uv, 2)
    return uv[1:2:n], uv[2:2:n]
end

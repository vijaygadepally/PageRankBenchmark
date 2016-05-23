#
# Read Edges from file
#
function StrFileRead(fname)
    uv = readdlm(fname, Int)
    return uv[1:2:end], uv[2:2:end]
end

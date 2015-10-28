#
# Read Edges from file
#
function StrFileRead(fname)
  fid = open(fname,"r");
  uv = split(readall(fid));
  # ut = int(uv[1:2:end])';   # makt it horizontal vector
  # vt = int(uv[2:2:end])';   # makt it horizontal vector
  ut = [parse(Int64,s) for s=uv[1:2:end]] # makt it horizontal vector
  vt = [parse(Int64,s) for s=uv[2:2:end]] # makt it horizontal vector
  close(fid)
  return ut,vt
end


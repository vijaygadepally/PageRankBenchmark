#
# Read Edges from file
#
function EdgeFileRead(fname)
  uv = readdlm(fname,'\t');
  ut = uv[:,1]';
  vt = uv[:,2]';
  return ut,vt
end


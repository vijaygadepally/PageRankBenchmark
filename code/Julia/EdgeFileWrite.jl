function EdgeFileWrite(u,v,fname)
#
# Write u and v pair per line to fname
#
  fid = open(fname,"w");
  writedlm(fid,[u; v].'); 
  close(fid);
end


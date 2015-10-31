function StrFileWrite(edgeStr,fname)
#
# Write edgeStr to fname
#
  fid = open(fname,"w");
  write(fid,edgeStr);    # space separated
  close(fid);
end


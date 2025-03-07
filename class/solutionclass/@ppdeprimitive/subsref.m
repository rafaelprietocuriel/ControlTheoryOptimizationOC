function out=subsref(ppdePrim,index)
%
%
try
    out=subsref(struct(ppdePrim),index);
    return
end
try
    out=subsref(struct(ppdePrim.ppdetrajectory),index);
catch
    lasterr;
end

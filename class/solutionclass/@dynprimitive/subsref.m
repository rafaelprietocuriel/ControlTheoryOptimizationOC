function out=subsref(dynPrim,index)
%
%
try
    out=subsref(struct(dynPrim),index);
    return
end
try
    out=subsref(struct(dynPrim.octrajectory),index);
catch
    lasterr;
end
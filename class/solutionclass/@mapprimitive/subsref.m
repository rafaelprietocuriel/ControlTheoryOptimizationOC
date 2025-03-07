function out=subsref(mapPrim,index)
%
%
try
    out=subsref(struct(mapPrim),index);
    return
end
try
    out=subsref(struct(mapPrim.doctrajectory),index);
catch
    lasterr;
end
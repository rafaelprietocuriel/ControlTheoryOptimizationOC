function out=subsref(gdynPrim,index)
%
%

try
    out=subsref(struct(gdynPrim),index);
catch
    out=subsref(struct(octrajectory(gdynPrim)),index);
end

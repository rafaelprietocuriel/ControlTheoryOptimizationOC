function out=subsref(ocAsym,index)
%
%

try
    out=subsref(struct(ocAsym),index);
catch
    out=subsref(struct(octrajectory(ocAsym)),index);
end

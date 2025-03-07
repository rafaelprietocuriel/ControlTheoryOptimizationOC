function out=subsref(docAsym,index)
%
%

try
    out=subsref(struct(docAsym),index);
catch
    out=subsref(struct(doctrajectory(docAsym)),index);
end

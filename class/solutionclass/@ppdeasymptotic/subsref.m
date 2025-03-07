function out=subsref(ppdeAsym,index)
%
%

try
    out=subsref(struct(ppdeAsym),index);
catch
    out=subsref(struct(ppdeAsym.ppdetrajectory),index);
end

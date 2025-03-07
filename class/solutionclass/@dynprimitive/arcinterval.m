function out=arcinterval(dynPrim,idx)

if nargin==1
    out=arcinterval(dynPrim.octrajectory);
else
    out=arcinterval(dynPrim.octrajectory,idx);
end


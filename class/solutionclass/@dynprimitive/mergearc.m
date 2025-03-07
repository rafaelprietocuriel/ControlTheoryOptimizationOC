function dynPrim=mergearc(dynPrim,idx,varargin)

if nargin==1
    idx=[];
end
if isperiodic(dynPrim)
    dynPrim.octrajectory=mergearc(dynPrim.octrajectory,idx,varargin{:});
end
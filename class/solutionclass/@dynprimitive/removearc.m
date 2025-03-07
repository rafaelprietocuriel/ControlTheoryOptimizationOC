function dynPrim=removearc(dynPrim,removearcidx,varargin)

if isperiodic(dynPrim)
    dynPrim.octrajectory=removearc(dynPrim.octrajectory,removearcidx,varargin{:});
end
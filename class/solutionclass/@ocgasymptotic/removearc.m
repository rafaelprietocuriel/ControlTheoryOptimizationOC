function ocgAsym=removearc(ocgAsym,removearcidx,varargin)

if isempty(ocgAsym) || isempty(removearcidx)
    return
end

ocgAsym.ocgtrajectory=removearc(ocgAsym.ocgtrajectory,removearcidx,varargin{:});

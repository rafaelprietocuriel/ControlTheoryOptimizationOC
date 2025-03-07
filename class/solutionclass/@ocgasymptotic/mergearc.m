function ocgAsym=mergearc(ocgAsym,idx,varargin)

if isempty(ocgAsym) || isempty(idx)
    return
end

ocgAsym.ocgtrajectory=mergearc(ocgAsym.ocgtrajectory,idx,varargin{:});

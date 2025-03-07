function ocAsymNew=removearc(ocAsym,removearcidx,varargin)

if isempty(ocAsym) || isempty(removearcidx)
    ocAsymNew=ocAsym;
    return
end
ocAsymNew=ocasymptotic(removearc(octrajectory(ocAsym),removearcidx,varargin{:}),limitset(ocAsym));
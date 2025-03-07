function ocAsymNew=mergearc(ocAsym,idx,varargin)

if isempty(ocAsym) || isempty(idx)
    ocAsymNew=ocAsym;
    return
end
ocAsymNew=ocasymptotic(mergearc(octrajectory(ocAsym),idx,varargin{:}),limitset(ocAsym));
function ocTrjNew=removearc(ocTrj,removearcidx,varargin)

if isempty(ocTrj) || isempty(removearcidx)
    ocTrjNew=ocTrj;
    return
end
structTrj=struct(ocTrj);
structTrj.solver='';
%structTrj.solverinfo=[];
structTrj.solverinfo.pathtype=pathtype(ocTrj);
ocTrjNew=octrajectory(removearc(structTrj,removearcidx,varargin{:}));
function ocTrjNew=changearc(ocTrj,changearcidx,newarcid,varargin)

if isempty(ocTrj) || isempty(changearcidx)
    ocTrjNew=ocTrj;
    return
end
structTrj=struct(ocTrj);
structTrj.solver='';
structTrj.solverinfo=[];
structTrj.solverinfo.pathtype=pathtype(ocTrj);
ocTrjNew=octrajectory(changearc(structTrj,changearcidx,newarcid,varargin{:}));
function ocTrjNew=mergearc(ocTrj,idx,varargin)

if isempty(ocTrj)
    ocTrjNew=ocTrj;
    return
end
structTrj=struct(ocTrj);
structTrj.solver='';
%structTrj.solverinfo=[];
structTrj.solverinfo.pathtype=pathtype(ocTrj);
ocTrjNew=octrajectory(mergearc(structTrj,idx,varargin{:}));
function hocTrj=mergearc(hocTrj,idx,varargin)

arcarg=[];
if nargin>=3
    arcarg=varargin{1};
end
if isempty(hocTrj) || isempty(idx)
    return
end
if ~isempty(arcarg)
    hocTrj.arcarg(idx)=arcarg;
end
structTrj=struct(hocTrj);
structTrj.jumparg(idx(2))=[];
structTrj.solver='';
structTrj.solverinfo=[];
structTrj.solverinfo.pathtype=pathtype(hocTrj);

Xbd=structTrj.y(:,[1 end]);
structTrj.y(:,[1 end])=[];
structTrj.x([1 end])=[];
structTrj.arcposition=structTrj.arcposition-1;

structTrj=mergearc(structTrj,idx(1));

structTrj.y=[Xbd(:,1) structTrj.y Xbd(:,2)];
structTrj.x=structTrj.x([1 1:end end]);
structTrj.arcposition=structTrj.arcposition+1;
hocTrj=hybridoctrajectory(structTrj);

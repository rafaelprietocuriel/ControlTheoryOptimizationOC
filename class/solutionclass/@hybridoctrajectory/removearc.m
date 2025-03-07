function hocTrj=removearc(hocTrj,removearcidx,varargin)

if isempty(hocTrj) || isempty(removearcidx)
    return
end
structTrj=struct(hocTrj);
Xbd=structTrj.y(:,[1 end]);
structTrj.y(:,[1 end])=[];
structTrj.x([1 end])=[];
structTrj.arcposition=structTrj.arcposition-1;
structTrj.solverinfo=[];
structTrj.solverinfo.pathtype=pathtype(hocTrj);
structTrj=removearc(structTrj,removearcidx,varargin{:});
structTrj.jumparg(removearcidx)=[];
structTrj.y=[Xbd(:,1) structTrj.y Xbd(:,2)];
structTrj.x=structTrj.x([1 1:end end]);
structTrj.arcposition=structTrj.arcposition+1;
hocTrj=hybridoctrajectory(structTrj);
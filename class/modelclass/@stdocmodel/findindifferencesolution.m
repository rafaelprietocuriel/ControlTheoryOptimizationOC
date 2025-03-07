function ocMul=findindifferencesolution(ocObj,idx1,idx2,ipt,varargin)
%
% FINDINDIFFERENCESOLUTION find approximate indifference solution.

contRes=contresult(ocObj,[idx1,idx2]);

idx=[0 0];
for ii=1:length(contRes)
    contSol=contRes{ii}.ContinuationSolution;
    ocEP=contRes{ii}.LimitSet;
    dist=zeros(1,length(contSol));
    for jj=1:length(contSol)
        ocEx=ocasymptotic(octrajectory(contSol(jj)),ocEP);
        x=state(ocObj,ocEx);
        if jj==1
            refvec=x(:,1)-ipt;
            refvec=refvec/norm(refvec);
        end
        dist(jj)=sign(sum(refvec.*(x(:,1)-ipt)))*norm(x(:,1)-ipt);
    end
    tmp=cont2idx(dist,0);
    if ~isempty(tmp)
        idx(ii)=tmp(1);
    end
end

if any(idx==0)
    ocMul=[];
    return
end
contSol=contRes{1}.ContinuationSolution;
ocMul{1}=ocasymptotic(octrajectory(contSol(idx(1))),contRes{1}.LimitSet);
contSol=contRes{2}.ContinuationSolution;
ocMul{2}=ocasymptotic(octrajectory(contSol(idx(2))),contRes{2}.LimitSet);
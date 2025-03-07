function [X,par,idx,arcarg]=matcontsolution2ep(matRes,searchdata,idx)
X=[];
par=[];
arcarg=[];
if isempty(matRes)
    return
end
if nargin==2
    idx=[];
end
if ischar(idx)
    if strcmp(idx,'end')
        idx=length(matRes.ContinuationSolution.userinfo.varyparametervalue);
    elseif strcmp(idx,'end-1')
        idx=length(matRes.ContinuationSolution.userinfo.varyparametervalue)-1;
    end
end
par0=matRes.ContinuationSolution.modelparameter;
par0=par0(:).';

if isempty(idx)
    switch matRes.ContinuationClassification
        case 'modelequilibrium'
            idx=cont2idx(matRes.ContinuationSolution.userinfo.varyparametervalue,searchdata);
        case 'modellimitpoint'
            idx=cont2idx(matRes.ContinuationSolution.userinfo.varyparametervalue(searchdata(1),:),searchdata(2));
    end
end

arcarg=zeros(1,length(idx));
par=repmat(par0,length(idx),1);
X=zeros(size(matRes.ContinuationSolution.y,1),length(idx));
for ii=1:length(idx)
    X(:,ii)=matRes.ContinuationSolution.y(:,idx(ii));
    par(ii,:)=par0;
    par(ii,matRes.ContinuationSolution.userinfo.varyparameterindex)=matRes.ContinuationSolution.userinfo.varyparametervalue(:,idx(ii));
    arcarg(ii)=matRes.ContinuationSolution.arcarg(idx(ii));
end

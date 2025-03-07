function [X,par,idx,arcarg]=matcontsolution2lp(matRes,searchdata,idx)
X=[];
par=[];
if isempty(matRes)
    return
end
if nargin==2
    idx=[];
end

par0=matRes.ContinuationSolution.modelparameter;

c=matRes.ContinuationInformation;
switch matRes.ContinuationClassification
    case 'modelequilibrium'
        ctr=0;
        idx=zeros(1,length(c));
        arcarg=zeros(1,length(c));
        par=repmat(par0,length(c),1);
        X=zeros(size(matRes.ContinuationSolution.y,1),length(c));
        for ii=1:length(c)
            if strcmp(c(ii).label,'LP')
                ctr=ctr+1;
                X(:,ctr)=matRes.ContinuationSolution.y(:,c(ii).index);
                par(ctr,:)=par0;
                par(ctr,matRes.ContinuationSolution.userinfo.varyparameterindex)=matRes.ContinuationSolution.userinfo.varyparametervalue(c(ii).index);
                arcarg(ctr)=matRes.ContinuationSolution.arcarg(c(ii).index);
                idx(ctr)=c(ii).index;
            end
        end
        X(:,ctr+1:end)=[];
        par(ctr+1:end,:)=[];
        arcarg(ctr+1:end)=[];
        idx(ctr+1:end)=[];
    case 'modellimitpoint'
        % searchdata; [parameterindex parametervalue]
        par0=matRes.ContinuationSolution.modelparameter(:).';
        if ischar(idx)
            if strcmp(idx,'end')
                idx=length(matRes.ContinuationSolution.userinfo.varyparametervalue);
            elseif strcmp(idx,'end-1')
                idx=length(matRes.ContinuationSolution.userinfo.varyparametervalue)-1;
            end
        elseif isempty(idx)
            if isempty(searchdata)
                idx=[];
                arcarg=[];
                return
            end
            idx=cont2idx(matRes.ContinuationSolution.userinfo.varyparametervalue(searchdata(1),:),searchdata(2));
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
end
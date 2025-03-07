function [X,s,par,arcid]=matcontsolution2bp(matRes,varargin)
X=[];
par=[];
s=[];

if isempty(matRes)
    return
end

par0=matRes.ContinuationSolution.modelparameter;
arcid=matRes.ContinuationSolution.arcarg(1);

c=matRes.ContinuationInformation;
switch matRes.ContinuationClassification
    case 'modelequilibrium'
        ctr=0;
        s=c;
        idx=[];
        for ii=1:length(c)
            if strcmp(c(ii).label,'BP')
                ctr=ctr+1;
                idx(ctr)=ii;
                X=[X matRes.ContinuationSolution.y(:,c(ii).index)];
                par(ctr,:)=par0;
                par(ctr,matRes.ContinuationSolution.userinfo.varyparameterindex)=matRes.ContinuationSolution.userinfo.varyparametervalue(c(ii).index);
            end
        end
        s=s(idx);
    case 'modellimitpoint'
        searchdata=[];
        if nargin>=2
            searchdata=varargin{1}; % [parameterindex parametervalue]
        end
        if isempty(searchdata)
            return
        end
        idx=cont2idx(matRes.ContinuationSolution.userinfo.varyparametervalue(searchdata(1),:),searchdata(2));
        
        for ii=1:length(idx)
                X=[X matRes.ContinuationSolution.y(:,idx(ii))];
                par(ii,:)=par0;
                par(ii,matRes.ContinuationSolution.userinfo.varyparameterindex)=matRes.ContinuationSolution.userinfo.varyparametervalue(:,idx(ii));
        end
end
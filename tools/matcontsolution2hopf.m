function [X,par]=matcontsolution2hopf(matRes,varargin)
X=[];
par=[];

if isempty(matRes)
    return
end

par0=matRes.ContinuationSolution.modelparameter;

c=matRes.ContinuationInformation;
switch matRes.ContinuationClassification
    case 'modelequilibrium'
        ctr=0;
        for ii=1:length(c)
            if strcmp(strtrim(c(ii).label),'H') && strcmp(strtrim(c(ii).msg),'Hopf')
                ctr=ctr+1;
                X=[X matRes.ContinuationSolution.y(:,c(ii).index)];
                par(ctr,:)=par0;
                par(ctr,matRes.ContinuationSolution.userinfo.varyparameterindex)=matRes.ContinuationSolution.userinfo.varyparametervalue(c(ii).index);
            end
        end
    case 'modelhopf'
        searchdata=[];
        if nargin>=2
            searchdata=varargin{1}; % [parameterindex parametervalue]
        end
        if isempty(searchdata)
            return
        end
        if length(searchdata)==1
            idx=searchdata;
        else
            idx=cont2idx(matRes.ContinuationSolution.userinfo.varyparametervalue(searchdata(1),:),searchdata(2));
        end
        for ii=1:length(idx)
            X=[X matRes.ContinuationSolution.y(:,idx(ii))];
            par(ii,:)=par0;
            par(ii,matRes.ContinuationSolution.userinfo.varyparameterindex)=matRes.ContinuationSolution.userinfo.varyparametervalue(:,idx(ii));
        end
end
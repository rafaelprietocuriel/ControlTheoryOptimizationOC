function [ocLC,par]=matcontsolution2lc(matRes,parval,paridx)
ocLC=[];
par=[];

if isempty(matRes)
    return
end

par0=matRes.ContinuationSolution(1).modelparameter;
par0=par0(:).';

switch matRes.ContinuationClassification
    case 'modellimitcycle'
        if nargin==2
            idx=parval;
        else
            actparidx=find(matRes.ContinuationSolution(end).modelparameter-matRes.ContinuationSolution(floor(length(matRes.ContinuationSolution)/2)).modelparameter);
            n=length(matRes.ContinuationSolution);
            par=zeros(1,n);
            for ii=1:n
                par(ii)=matRes.ContinuationSolution(ii).modelparameter(actparidx);
            end
            idx=cont2idx(par,parval);
        end
    case 'modellimitpointcycle'
        if nargin==2
            idx=parval;
        else
            actparidx=find(matRes.ContinuationSolution(end).modelparameter-matRes.ContinuationSolution(1).modelparameter);
            actparidx0=actparidx(paridx);
            n=length(matRes.ContinuationSolution);
            par=zeros(1,n);
            for ii=1:n
                par(ii)=matRes.ContinuationSolution(ii).modelparameter(actparidx0);
            end
            idx=cont2idx(par,parval);
        end
end

par=[];
ctr=0;
par=zeros(length(idx),length(par0));
for ii=1:length(idx)
    ctr=ctr+1;
    ocLCStruct.octrajectory=octrajectory(matRes.ContinuationSolution(idx(ii)));
    ocLCStruct.period=matRes.ContinuationSolution(ii).timehorizon;
    if length(idx)>1
        ocLC{ctr}=dynprimitive(ocLCStruct);
    else
        ocLC=dynprimitive(ocLCStruct);
    end
    %par(ctr,:)=par0;
    par(ctr,:)=matRes.ContinuationSolution(idx(ii)).modelparameter;
end

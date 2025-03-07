function [idx minval]=cont2idx(ocObj,contindex,cmpcoord,cmpval)
%
% CONT2IDX returns indices where a curve crosses a specific parameter value
%
% IDX=CONT2IDX(CONTPAR,PARVAL) returns the index/indices where the
% curve stored as a vector of discrete numbers CONTPAR, usually coming from
% a continuation process, crosses the parameter value PARVAL.

contresultStruct=contresult(ocObj);
numresult=numel(contresultStruct);
if contindex<1 || contindex>max(numresult)
    ocmaterror('')
end
stepnum=numel(contresultStruct{contindex}.ContinuationSolution);
yt=zeros(numel(cmpcoord),stepnum);
for ii=1:stepnum
    yt(:,ii)=contresultStruct{contindex}.ContinuationSolution(ii).y(cmpcoord,1);
end

[minval idx]=min(sqrt(sum((yt-cmpval(:,ones(1,stepnum))).^2)));
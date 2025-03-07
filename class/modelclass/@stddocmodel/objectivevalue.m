function o=objectivevalue(docObj,solObj,varargin)
%
% 
o=[];
if isempty(docObj)
    return
end
if nargin==1
    return
end

o=sum(objectivefunction(docObj,solObj,1));
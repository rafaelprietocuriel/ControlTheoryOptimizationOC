function [ocEx idx]=nonadmissiblesolution(ocObj,idx)
%
% EXTREMALSOLUTION returns the slice manifold(s) stored in OCOBJ
%
% OCEX=EXTREMALSOLUTION(OCOBJ) returns a cell arry of the extremal solutions
% stored in the  'Result' field 'Continuation.ExtremalSolution' of OCOBJ.

% variable declaration
ocEx=occurve();

contResultStruct=contresult(ocObj);

if isempty(ocObj) || isempty(contResultStruct)
    return
end

if nargin==1 || isempty(idx)
    idx=1:numel(contResultStruct);
end
f=cellfun(@isfield,contResultStruct(idx),repmat({'NonadmissibleSolution'},1,numel(idx)),'UniformOutput',false);
idx(~[f{:}])=[];
ocEx=cellfun(@getfield,contResultStruct(idx),repmat({'NonadmissibleSolution'},1,numel(idx)),'UniformOutput',false);
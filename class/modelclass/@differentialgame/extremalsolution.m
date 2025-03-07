function [ocEx idx]=extremalsolution(dgObj,idx,varargin)
%
% EXTREMALSOLUTION returns the slice manifold(s) stored in OCOBJ
%
% OCEX=EXTREMALSOLUTION(OCOBJ) returns a cell arry of the extremal solutions
% stored in the  'Result' field 'Continuation.ExtremalSolution' of OCOBJ.

% variable declaration
ocEx=occurve();
contname='';

if nargin>2
    contname=varargin{1};
end
if isempty(contname)
    contResultStruct=contresult(dgObj);
else
    contResultStruct=result(ocObj,contname);
end

if isempty(dgObj) || isempty(contResultStruct)
    return
end

if nargin==1 || isempty(idx)
    idx=1:numel(contResultStruct);
end
f=cellfun(@isfield,contResultStruct(idx),repmat({'ExtremalSolution'},1,numel(idx)),'UniformOutput',false);
idx(~[f{:}])=[];
ocEx=cellfun(@getfield,contResultStruct(idx),repmat({'ExtremalSolution'},1,numel(idx)),'UniformOutput',false);
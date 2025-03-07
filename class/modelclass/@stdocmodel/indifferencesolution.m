function [ocEx idx]=indifferencesolution(ocObj,idx)
%

% variable declaration
ocEx=occurve();

contResultStruct=contresult(ocObj);

if isempty(ocObj) || isempty(contResultStruct)
    return
end

if nargin==1 || isempty(idx)
    idx=1:numel(contResultStruct);
end
f=cellfun(@isfield,contResultStruct(idx),repmat({'IndifferencePointCurve'},1,numel(idx)),'UniformOutput',false);
idx(~[f{:}])=[];
ocEx=cellfun(@getfield,contResultStruct(idx),repmat({'IndifferencePointCurve'},1,numel(idx)),'UniformOutput',false);
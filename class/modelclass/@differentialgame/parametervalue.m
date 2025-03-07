function [parval,parvar]=parametervalue(dgObj,varargin)
%
% PARAMETERVALUE returns the set of parameter values for the actual
% oc-model.
%
% [PARVAL,PARVAR]=PARAMETERVALUE(OCOBJ) OCOBJ is an instance of the
% stdocmodel class. The function returns the values of the parameter PARVAL
% and a cell array of strings consisting of the parameter variable names
% PARVAR.

parval=[];
parvar=[];
parindex=[];
if isempty(dgObj)
    return
end
if nargin==2
    parindex=parameterindex(dgObj,varargin{1});
end
info=retrievedifferentialgameinformation(dgObj.Model,'parametervalue');
parval=info.value;
if nargout==2
    info=retrievedifferentialgameinformation(dgObj.Model,'parametername');
    parvar=[info.value(:).'];
end

if ~isempty(parindex)
    parval=parval(parindex);
    if ~isempty(parvar)
        parvar=parvar(parindex);
    end
end
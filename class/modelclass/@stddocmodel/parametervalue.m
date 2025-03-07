function [parval,parvar]=parametervalue(ocObj,varargin)
%
% PARAMETERVALUE returns the set of parameter values for the actual
% discrete ocmodel.
%
% [PARVAL,PARVAR]=PARAMETERVALUE(OCOBJ) OCOBJ is an instance of the
% stddocmodel class. The function returns the values of the parameter PARVAL
% and a cell array of strings consisting of the parameter variable names
% PARVAR.

parval=[];
parvar=[];
parindex=[];
if isempty(ocObj)
    return
end
if nargin==2
    parindex=parameterindex(ocObj,varargin{1});
end
info=retrievediffmodelinformation(ocObj.Model,'parametervalue');
parval=info.value;
if nargout==2
    info=retrievediffmodelinformation(ocObj.Model,'parametername');
    parvar=[info.value(:).'];
end

if ~isempty(parindex)
    parval=parval(parindex);
    if ~isempty(parvar)
        parvar=parvar(parindex);
    end
end
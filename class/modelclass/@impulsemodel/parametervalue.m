function [parval,parvar]=parametervalue(ocObj,varargin)
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
if isempty(ocObj)
    return
end
if nargin==2
    if ~isempty(varargin{1})
        parindex=parameterindex(ocObj,varargin{1});
        if isempty(parindex)
            return
        end
    else
        parindex=parameterindex(ocObj);
    end
end
info=retrievemodelinformation(ocObj,'parametervalue');
parval=info.value;
if nargout==2
    info=retrievemodelinformation(ocObj,'parametername');
    parvar=[info.value(:).'];
end

if ~isempty(parindex)
    parval=parval(parindex);
    if ~isempty(parvar)
        parvar=parvar(parindex);
    end
end
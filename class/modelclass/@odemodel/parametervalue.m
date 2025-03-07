function [parval,parvar]=parametervalue(odeObj,varargin)
%
% PARAMETERVALUE returns the set of parameter values for the actual
% ode-model.
%
% [PARVAL,PARVAR]=PARAMETERVALUE(ODEOBJ) ODEOBJ is an instance of the
% stdocmodel class. The function returns the values of the parameter PARVAL
% and a cell array of strings consisting of the parameter variable names
% PARVAR.

parval=[];
parvar=[];
parindex=[];
if isempty(odeObj)
    return
end
if nargin==2
    parindex=parameterindex(odeObj,varargin{1});
end
info=retrievemodelinformation(odeObj.Model,'parametervalue');
parval=info.value;
if nargout==2
    info=retrievemodelinformation(odeObj.Model,'parametername');
    parvar=[info.value(:).'];
end

if ~isempty(parindex)
    parval=parval(parindex);
    if ~isempty(parvar)
        parvar=parvar(parindex);
    end
end
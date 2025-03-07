function [parval,parvar]=parametervalue(ppdeocObj,varargin)
%
% PARAMETERVALUE returns the set of parameter values for the actual
% ppdeoc-model.
%
% [PARVAL,PARVAR]=PARAMETERVALUE(PPDEOCOBJ) PPDEOCOBJ is an instance of the
% stdocmodel class. The function returns the values of the parameter PARVAL
% and a cell array of strings consisting of the parameter variable names
% PARVAR.

parval=[];
parvar=[];
parindex=[];
if isempty(ppdeocObj)
    return
end
if nargin==2
    parindex=parameterindex(ppdeocObj,varargin{1});
end
info=retrieveppdemodelinformation(ppdeocObj.Model,'parametervalue');
parval=info.value;
if nargout==2
    info=retrieveppdemodelinformation(ppdeocObj.Model,'parametername');
    parvar=info.value(:).';
end

if ~isempty(parindex)
    parval=parval(parindex);
    if ~isempty(parvar)
        parvar=parvar(parindex);
    end
end
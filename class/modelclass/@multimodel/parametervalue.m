function [parval,parvar]=parametervalue(mmObj,varargin)
%
% PARAMETERVALUE returns the set of parameter values for the actual
% oc-model.
%
% [PARVAL,PARVAR]=PARAMETERVALUE(OCOBJ) OCOBJ is an instance of the
% stdocmodel class. The function returns the values of the parameter PARVAL
% and a cell array of strings consisting of the parameter variable names
% PARVAR.

nummod=numberofmodels(mmObj);
parval=cell(1,nummod);
parvar=cell(1,nummod);
if isempty(mmObj)
    return
end
if nargin==1
    inparvar=[];
end
if nargin==2
    inparvar=varargin{1};
end
for ii=1:nummod
    if iscell(inparvar)
        if isempty(inparvar{ii})
            [parval{ii},parvar{ii}]=parametervalue(mmObj.Model{ii});
        else
            [parval{ii},parvar{ii}]=parametervalue(mmObj.Model{ii},inparvar{ii});
        end
    else
        if isempty(inparvar)
            [parval{ii},parvar{ii}]=parametervalue(mmObj.Model{ii});
        else
            [parval{ii},parvar{ii}]=parametervalue(mmObj.Model{ii},inparvar);
        end
    end
end

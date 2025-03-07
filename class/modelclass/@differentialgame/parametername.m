function parvar=parametername(dgObj,varargin)
%
% PARAMETERNAME returns the set of parameter values for the actual
% oc-model.
%
% [PARVAL,PARVAR]=PARAMETERVALUE(dgObj) dgObj is an instance of the
% stdocmodel class. The function returns the values of the parameter PARVAL
% and a cell array of strings consisting of the parameter variable names
% PARVAR.

parvar=[];
parindex=[];
if isempty(dgObj)
    return
end
if nargin==2
    if ~isempty(varargin{1})
        parindex=varargin{1};
        if isempty(parindex)
            return
        end
    else
        return
    end
end
    
info=retrievemodelinformation(dgObj.Model,'parametername');
parvar=info.value;
if nargin==1
    return
end
parvar=parvar(parindex);
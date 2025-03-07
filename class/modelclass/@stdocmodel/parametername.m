function parvar=parametername(ocObj,varargin)
%
% PARAMETERNAME returns the set of parameter values for the actual
% oc-model.
%
% [PARVAL,PARVAR]=PARAMETERVALUE(OCOBJ) OCOBJ is an instance of the
% stdocmodel class. The function returns the values of the parameter PARVAL
% and a cell array of strings consisting of the parameter variable names
% PARVAR.

parvar=[];
parindex=[];
if isempty(ocObj)
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
    
info=retrievemodelinformation(ocObj.Model,'parametername');
parvar=info.value;
if nargin==1
    return
end
parvar=parvar(parindex);
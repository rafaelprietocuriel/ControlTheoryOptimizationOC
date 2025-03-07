function out=salvagevalue(ocObj,solObj,varargin)
%
% SALVAGEVALUE

out=[];
if isempty(ocObj)
    out=[];
    return
end
if nargin==1
    solObj=[];
end
par=parametervalue(ocObj);
if ishybridoctrajectory(solObj)
    jumparg=jumpargument(solObj);
    arcarg=arcargument(solObj);
    indepvar=time(ocObj,solObj,1);
    depvar=dependentvar(solObj);
    jumparg=jumparg(end);
    arcarg=arcarg(end);
    indepvar=indepvar(end);
    depvar=depvar(:,end-1:end);
elseif isnumeric(solObj)
    arcarg=0;
    jumparg=0;
    indepvar=0;
    depvar=[solObj solObj];
end

if isempty(solObj)
    try
        out=feval(ocObj,'SymbolicSalvagevalue',arcarg);
    catch
        return
    end
else
    % return value of the canonical system evaluated at 'depvar'
    out=feval(ocObj,'Salvagevalue',indepvar,depvar,par,arcarg,jumparg);
end

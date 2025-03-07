function out=salvagevalue(ocObj,solObj,ct,varargin)
%
% EXOGENOUSFUNCTION returns the value of the term specified in the model
% specific UserFunction file.

if isempty(ocObj)
    out=[];
    return
end
par=parametervalue(ocObj);
if isoctrajectory(solObj) || isocasymptotic(solObj)
    arcarg=arcargument(solObj);
    indepvar=time(ocObj,solObj,1);
    depvar=dependentvar(solObj);
end

out=feval(ocObj,'Salvagevalue',indepvar(end),depvar(:,end),par,arcarg(end),ct);

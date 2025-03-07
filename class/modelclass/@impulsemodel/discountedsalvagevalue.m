function out=discountedsalvagevalue(ocObj,solObj,varargin)
%
% SALVAGEVALUE

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
    indepvar=time(solObj,1);
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
        varargout{1}=feval(ocObj,'DiscountedSalvagevalue',[],[],par,[]);
    catch
        varargout{1}=[];
        return
    end
else
    % return value of the canonical system evaluated at 'depvar'
    out=feval(ocObj,'DiscountedSalvagevalue',indepvar,depvar,par,arcarg,jumparg);
end

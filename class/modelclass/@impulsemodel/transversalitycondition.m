function out=transversalitycondition(ocObj,solObj,varargin)
%
% TRANSVERSALITYCONDITION

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
    endtime=arcinterval(ocObj,solObj,1);
    endtime=endtime(end);
    depvar=dependentvar(solObj);
    jumparg=jumparg(end);
    arcarg=arcarg(end);
    depvar=depvar(:,end-1:end);
elseif isstruct(solObj)
    arcarg=solObj.arcarg(end);
    jumparg=solObj.jumparg(end);
    endtime=solObj.arcinterval;
    endtime=endtime(end);
    depvar=solObj.y(:,end-1:end);
elseif isnumeric(solObj)
    arcarg=0;
    jumparg=0;
    endtime=0;
    depvar=[solObj solObj];
end
if isempty(solObj)
    try
        varargout{1}=feval(ocObj,'TransversalityBC',[],[],par,[]);
    catch
        varargout{1}=[];
        return
    end
else
    % return value of the canonical system evaluated at 'depvar'
    out=feval(ocObj,'TransversalityBC',endtime,depvar,par,arcarg,jumparg);
end

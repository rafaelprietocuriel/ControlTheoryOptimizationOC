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
if isnumeric(solObj)
    arcarg=0;
    indepvar=0;
    depvar=[solObj solObj];
elseif isdynprimitive(solObj) || isoctrajectory(solObj) || isocasymptotic(solObj)
    arcarg=arcargument(solObj);
    indepvar=time(ocObj,solObj,1);
    depvar=dependentvar(solObj);
elseif isoccurve(solObj)
    arcarg=arcargument(solObj);
    depvar=dependentvar(solObj);
    indepvar=zeros(1,size(depvar,2));
elseif isocgradtrajectory(solObj)
    switch octype(solOb)
        case 'concentrated'
            arcarg=[];
            indepvar=time(ocObj,solObj);
            depvar=dependentvar(solObj);
    end
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
    out=feval(ocObj,'TransversalityBC',indepvar(end),depvar(:,end),par,arcarg);
end

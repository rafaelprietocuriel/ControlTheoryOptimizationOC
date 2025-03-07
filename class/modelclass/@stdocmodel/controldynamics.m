function varargout=controldynamics(ocObj,solObj)
%
% DERIVATIVEGENERALIZEDCONTROL returns the derivative of the generalized
% control(s) (control(s) and Lagrangemultiplier with respect to the
% state(s) and costate(s).
%

if isempty(ocObj)
    ocmatmsg('Input model is empty.\n')
    varargout{1}=[];
    return
end
if ~isimplicit(ocObj)
    ocmatmsg('Input model has no implicit control.\n')
    varargout{1}=[];
    return
end
if nargin==1
    solObj=[];
end

if isempty(solObj)
    ocmatmsg('Input solution is empty.\n')
    varargout{1}=[];
    return
end

par=parametervalue(ocObj);
if isgdynprimitive(solObj)
    arcarg=arcargument(solObj);
    arcn=arcnum(solObj);
    indepvar=time(ocObj,solObj,1);
    depvar=dependentvar(solObj);
else
    ocmatmsg('Not implemented for class %s.\n',class(solObj));
    return
end
dudt=cell(1,arcn);
for ii=1:arcn
    if ~isimplicit(ocObj,arcarg(ii))
        ocmatmsg('Solution object has no implicit control.\n')
        varargout{1}=[];
        return
    end
    dudt{ii}=feval(ocObj,'OptimalControlDynamics',indepvar,depvar{ii},par,arcarg(ii));
end
varargout{1}=dudt;

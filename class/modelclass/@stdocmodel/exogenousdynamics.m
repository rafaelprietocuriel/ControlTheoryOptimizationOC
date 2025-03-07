function varargout=exogenousdynamics(ocObj,solObj,varargin)
%
% EXOGENOUSFUNCTION returns the value of the term specified in the model
% specific UserFunction file.

if isempty(ocObj)
    varargout{1}=[];
    return
end
par=parametervalue(ocObj);
if isoctrajectory(solObj) || isocasymptotic(solObj)
    arcarg=arcargument(solObj);
    indepvar=time(ocObj,solObj,1);
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    arcn=arcnum(solObj);
end

% return value of the canonical system evaluated at 'depvar'
connectflag=[];
if nargin>=3
    connectflag=varargin{1};
end
if isempty(connectflag)
    connectflag=0;
end
for ii=1:arcn
    arcp=arcpos(1,ii):arcpos(2,ii);
    out=feval(ocObj,'ExogenousDynamics',indepvar(arcp),depvar(:,arcp),par,arcarg(ii));
    if connectflag
        varargout{1}(:,arcp)=out;
    else
        varargout{ii}=out;
    end
end

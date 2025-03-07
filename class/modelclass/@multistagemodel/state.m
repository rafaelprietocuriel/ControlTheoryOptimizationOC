function varargout=state(ocObj,solObj,ct,varargin)
%
% USERFUNCTION returns the value of the term specified in the model
% specific UserFunction file.

if isempty(ocObj)
    varargout{1}=[];
    return
end
par=parametervalue(ocObj);
if isoctrajectory(solObj) || isocasymptotic(solObj)
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    arcn=arcnum(solObj);
end

% return value of the canonical system evaluated at 'depvar'
connectflag=[];
if nargin>=4
    connectflag=varargin{1};
end
if isempty(connectflag)
    connectflag=0;
end
for ii=1:arcn
    arcp=arcpos(1,ii):arcpos(2,ii);
    xcoord=statecoord(ocObj);
    lcoord=statecoord(ocObj);
    l=depvar(lcoord,arcp);
    if connectflag
        varargout{1}(xcoord,arcp)=l;
    else
        varargout{ii}=l;
    end
end

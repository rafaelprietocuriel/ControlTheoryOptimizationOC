function varargout=time(ocObj,solObj,varargin)
%

connectflag=[];
if isempty(ocObj)
    return
end
if nargin>=3
    connectflag=varargin{1};
end
if isempty(connectflag)
    connectflag=0;
end
if ishybridoctrajectory(solObj)
    indepvar=independentvar(solObj);
    arcpos=arcposition(solObj);
    arcint=arcinterval(solObj);
    x0=initialtime(solObj);
    arcn=arcnum(solObj);
else
    return
end

diffarcinterval=diff(arcint);
t=x0;
for ii=1:arcn
    arcp=arcpos(1,ii):arcpos(2,ii);
    t=t(end)+(indepvar(1,arcp)-ii+1)*diffarcinterval(ii);
    if connectflag
        varargout{1}(1,arcp)=t;
    else
        varargout{ii}=t;
    end
end
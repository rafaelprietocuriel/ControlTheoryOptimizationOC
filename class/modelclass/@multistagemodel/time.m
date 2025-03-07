function out=time(mmObj,solObj,varargin)
%
% 

out=[];
connectflag=[];
if isempty(mmObj)
    return
end
stageindex=stage(mmObj);
if nargin>=3
    connectflag=varargin{1};
end
if isempty(connectflag)
    connectflag=0;
end

if ismmultipath(solObj)
    if stageindex>numberofparts(solObj)
        return
    end
    indepvar=independentvar(solObj(stageindex));
    arcpos=arcposition(solObj(stageindex));
    arcint=arcinterval(solObj(stageindex));
    arcn=arcnum(solObj(stageindex));
elseif isoctrajectory(solObj)
    indepvar=independentvar(solObj);
    arcpos=arcposition(solObj);
    arcint=arcinterval(solObj);
    arcn=arcnum(solObj);
else
    return
end

% return optimal control value evaluated at 'depvar'
diffarcinterval=diff(arcint);

t=arcint(1);
for ii=1:arcn
    arcp=arcpos(1,ii):arcpos(2,ii);
    t=t(end)+(indepvar(1,arcp)-ii+1)*diffarcinterval(ii);
    if connectflag
        out(arcp)=t;
    else
        out{ii}=t;
    end
end

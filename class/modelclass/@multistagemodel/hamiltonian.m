function out=hamiltonian(mmObj,solObj,varargin)
%
% 

out=[];
if isempty(mmObj)
    return
end
par=parametervalue(mmObj);
stageindex=stage(mmObj);
if ~ismmultipath(solObj)
    return
end

if stageindex>numberofparts(solObj)
    return
end
arcarg=arcargument(solObj(stageindex));
t=time(mmObj,solObj(stageindex),1);
depvar=dependentvar(solObj(stageindex));
arcpos=arcposition(solObj(stageindex));
arcn=arcnum(solObj(stageindex));
ct=connectiontime(solObj);
% return optimal control value evaluated at 'depvar'
connectflag=[];
if nargin>=3
    connectflag=varargin{1};
end
if isempty(connectflag)
    connectflag=0;
end
for ii=1:arcn
    arcp=arcpos(1,ii):arcpos(2,ii);
    o=feval(mmObj,'Hamiltonian',t(arcp),depvar(:,arcp),par,arcarg(ii),ct);
    if connectflag
        out(arcp)=o;
    else
        out{ii}=o;
    end
end

function o=integralobjectivevalue(ocObj,solObj,varargin)
%
%
o=[];
indepvar=0;
if isempty(ocObj) || isempty(solObj)
    return
end
if nargin<=2
    arcarg=0;
end
if nargin==1
    solObj=[];
end
par=parametervalue(ocObj);
if ishybridoctrajectory(solObj)
    arcarg=arcargument(solObj);
    indepvar=time(solObj,1);
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    arcn=arcnum(solObj);
elseif isstruct(solObj)
    arcarg=solObj.arcarg;
    indepvar=solObj.x*solObj.arcinterval(end);
    depvar=solObj.y;
    if isfield(solObj,'arcposition')
        arcpos=solObj.arcposition;
    else
        arcpos=[1;length(indepvar)];
    end
    arcn=length(arcarg);
end

if size(depvar,1)==2*statenum(ocObj)+1
    o=depvar(end,end);
else
    of=objectivefunction(ocObj,solObj,1);
    o=sum(diff(indepvar).*(of(1:end-1)+of(2:end))/2);
%     for ii=1:arcn
%         arcp=arcpos(1,ii):arcpos(2,ii);
%         
%         %of=feval(ocObj,'ObjectiveFunction',indepvar(arcp-1),depvar(:,arcp),par,arcarg(ii),[]);
%         o=o+sum(diff(indepvar).*(of(1:end-1)+of(2:end))/2);
%         %o(arcp)=[o(end) o(end)+cumsum(diff(indepvar(arcp-1)).*(of(1:end-1)+of(2:end))/2)];
%     end
end
function out=objectivefunction(mmObj,solObj,varargin)
%
% 

out=[];
connectflag=[];
if isempty(mmObj)
    return
end
if ~ismmultipath(solObj)
    return
end
if nargin>=3
    connectflag=varargin{1};
end
if isempty(connectflag)
    connectflag=0;
end

par=parametervalue(mmObj);

np=numberofparts(solObj);

if np>numberofmodels(mmObj)
    return
end
arcarg=arcargument(solObj);
indepvar=time(mmObj,solObj);
depvar=dependentvar(solObj);
ct=connectiontime(solObj);
for part=1:numberofparts(solObj)
    arcn=arcnum(solObj(part));
    arcpos=arcposition(solObj(part));
    o=[];
    for ii=1:arcn
        arcp=arcpos(1,ii):arcpos(2,ii);
        val=feval(mmObj.Model{part},'ObjectiveFunction',indepvar{part}(arcp),depvar{part}(:,arcp),par{part},arcarg{part}(ii),ct);
        o=[o val];
    end
    if ~connectflag
        out{ii}=o;
    else
        out=[out o];
    end
end

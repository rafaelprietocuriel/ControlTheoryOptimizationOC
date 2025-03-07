function varargout=time(docObj,solObj,varargin)
%
connectflag=[];
if isempty(docObj)
    return
end
if nargin==1
    solObj=[];
end
if nargin>=3
    connectflag=varargin{1};
end
if isempty(connectflag)
    connectflag=0;
end
if isdoctrajectory(solObj) %|| isocasymptotic(solObj)
    indepvar=[initialtime(solObj) independentvar(solObj)];
    arcpos=arcposition(solObj);
    arcn=arcnum(solObj);
else
    ocmaterror('Not defined for class ''%s''!\n',class(solObj))
end

for ii=1:arcn
    arcp=arcpos(1,ii):arcpos(2,ii)+1;
    if connectflag
        varargout{1}(xcoord,arcp(1:end-1))=indepvar(arcp(2:end));
    else
        varargout{ii}=indepvar(arcp(2:end));
    end
end

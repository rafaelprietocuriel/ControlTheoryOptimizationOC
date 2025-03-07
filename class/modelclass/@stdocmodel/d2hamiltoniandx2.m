function varargout=d2hamiltoniandx2(ocObj,solObj,varargin)
%


if isempty(ocObj)
    return
end
if isempty(solObj)
    return
end
par=parametervalue(ocObj);
arcarg=arcargument(solObj);
indepvar=time(ocObj,solObj,1);
depvar=dependentvar(solObj);
arcpos=arcposition(solObj);
arcn=arcnum(solObj);
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
    equationnum=canonicalsystemdimension(ocObj,arcarg(ii));
    for jj=1:size(arcp,2)
        H(:,:,jj)=feval(ocObj,'D2HamiltonianDX2',indepvar(arcp(:,jj)),depvar(:,arcp(:,jj)),par,arcarg(ii));
    end
    if connectflag
        varargout{1}(1:equationnum,arcp)=H;
    else
        varargout{ii}=H;
    end
end

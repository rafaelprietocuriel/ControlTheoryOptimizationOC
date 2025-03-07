function varargout=userfunction(ocObj,solObj,varargin)
%
% USERFUNCTION returns the value of the term specified in the model
% specific UserFunction file.

if isempty(ocObj)
    varargout{1}=[];
    return
end
if nargin==1
    solObj=[];
end
par=parametervalue(ocObj);
if isstruct(solObj)
    try
        arcarg=solObj.arcarg;
        indepvar=solObj.independentvar;
        depvar=solObj.dependentvar;
    catch
        ocmaterror('If the second input argument is a structure the fields ''dependentvar'', ''independentvar'' and ''arcarg'' have to exist!')
    end
elseif ishybridoctrajectory(solObj)
    arcarg=arcargument(solObj);
    jumparg=jumpargument(solObj);
    indepvar=time(ocObj,solObj,1);
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    arcn=arcnum(solObj);
elseif isoccurve(solObj)
    arcarg=arcargument(solObj);
    indepvar=zeros(size(arcarg));
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    if isempty(arcpos)
        arcarg=unique(arcarg);
        arcpos=[1 size(depvar,2)]';
        arcn=1;
    else
        arcpos=arcposition(solObj);
        arcn=1;
    end
elseif isnumeric(solObj)
    arcarg=0;
    indepvar=0;
    depvar=solObj;
    arcpos=[1;size(depvar,2)];
    arcn=1;
end

if isempty(solObj)
    try
        varargout{1}=feval(ocObj,'UserFunction',[],[],par,[]);
    catch
        varargout{1}=[];
        return
    end
else
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
        u=feval(ocObj,'UserFunction',indepvar(arcp-1),depvar(:,arcp),par,arcarg(ii),jumparg(ii));
        if connectflag
            varargout{1}(:,arcp-1)=u;
        else
            varargout{ii}=u;
        end
    end
end

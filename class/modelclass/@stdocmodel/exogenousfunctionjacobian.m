function varargout=exogenousfunctionjacobian(ocObj,solObj,exfuncnum,varargin)
%
% 

indepvar=0;
if isempty(ocObj)
    return
end
if nargin<=3
    arcarg=0;
end
if nargin==1
    solObj=[];
end
if isstruct(solObj)
    try
        arcarg=solObj.arcarg;
        try
            indepvar=solObj.independentvar;
        catch
            indepvar=solObj.x*solObj.arcinterval(end);
        end
        try
            depvar=solObj.dependentvar;
        catch
            depvar=solObj.y;
        end
        if isfield(solObj,'arcposition')
            arcpos=solObj.arcposition;
        else
            arcpos=[1;length(indepvar)];
        end
        arcn=numel(arcarg);
    catch
        ocmaterror(['Field is missing.\n' lasterr])
    end
elseif isdynprimitive(solObj) || isoctrajectory(solObj) || isocasymptotic(solObj)
    arcarg=arcargument(solObj);
    indepvar=time(ocObj,solObj,1);
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    arcn=arcnum(solObj);
elseif  isoccurve(solObj)
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
        arcn=arcnum(solObj);
    end
end
% return optimal control value evaluated at 'depvar'
connectflag=[];
if nargin>=1
    connectflag=varargin{1};
end
if isempty(connectflag)
    connectflag=0;
end
for ii=1:arcn
    arcp=arcpos(1,ii):arcpos(2,ii);
    equationnum=canonicalsystemdimension(ocObj,arcarg(ii));
    J=zeros(equationnum*exfuncnum,size(arcp,2));
    for jj=1:size(arcp,2)
        tmp=feval(ocObj,'ExogenousJacobian',indepvar(arcp(jj)),depvar(:,arcp(jj)),parametervalue(ocObj),arcarg(ii));
        if size(arcp,2)>1
            J(:,jj)=tmp(:);
        else
            J=tmp;
        end
    end
    if connectflag
        varargout{1}(1:equationnum*exfuncnum,arcp)=J;
    else
        varargout{ii}=J;
    end
end


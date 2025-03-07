function varargout=objectivefunction(dgObj,solObj,varargin)
%
% 

indepvar=0;
daeflag=0;
if isempty(dgObj)
    return
end
if nargin<=2
    arcarg=0;
end
if nargin==1
    solObj=[];
end
par=parametervalue(dgObj);
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
elseif isdaeoctrajectory(solObj)
    arcarg=[];
    indepvar=independentvar(solObj);
    depvar=dependentvar(solObj);
    daeflag=1;
elseif isdynprimitive(solObj) || isoctrajectory(solObj) || isocasymptotic(solObj)
    arcarg=arcargument(solObj);
    indepvar=time(dgObj,solObj,1);
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
if isempty(solObj)
    % return symbolic expressions
    % if repflag=1 control values/Lagrange multipliers are replaced by the
    % corresponding symbolic expressions for the optimal values
    repflag=[];
    arcarg=[];
    if nargin>=3
        arcarg=varargin{1};
    end
    if nargin>=4
        repflag=varargin{2};
    end
    if isempty(arcarg)
        arcarg=0;
    end
    if isempty(repflag)
        repflag=0;
    end
    varargout{1}=feval(dgObj,'SymbolicObjectiveFunction',arcarg,repflag);
else
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
        o=feval(dgObj,'ObjectiveFunction',indepvar(arcp),depvar(:,arcp),par,arcarg(ii));
        if connectflag
            varargout{1}(:,arcp)=o;
        else
            varargout{ii}=o;
        end
    end
end
function varargout=objectivefunction(ppdeObj,solObj,varargin)
%
% 

indepvar=0;
if isempty(ppdeObj)
    return
end
if nargin<=2
    arcarg=0;
end
if nargin==1
    solObj=[];
end
par=parametervalue(ppdeObj);
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
elseif ispdeprimitive(solObj) || ispdetrajectory(solObj) || ispdeasymptotic(solObj)
    arcarg=arcargument(solObj);
    indepvar=time(ppdeObj,solObj,1);
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    arcn=arcnum(solObj);
    femdat=femdata(solObj);
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
    varargout{1}=feval(ppdeObj,'SymbolicObjectiveFunction',arcarg,repflag);
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
        o=feval(ppdeObj,'ObjectiveFunction',indepvar(arcp),depvar(:,arcp),par,arcarg(ii),femdat);
        if connectflag
            varargout{1}(:,arcp)=o;
        else
            varargout{ii}=o;
        end
    end
end
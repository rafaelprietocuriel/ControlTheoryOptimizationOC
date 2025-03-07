function varargout=objectivefunction(docObj,solObj,varargin)
%
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
par=parametervalue(docObj);
if ismapprimitive(solObj) || isdoctrajectory(solObj) || isdocasymptotic(solObj)
    arcarg=arcargument(solObj);
    indepvar=[initialtime(solObj) independentvar(solObj)];
    depvar=[initialstate(solObj) dependentvar(solObj)];
    arcpos=arcposition(solObj);
    arcn=arcnum(solObj);
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
    varargout{1}=feval(docObj,'SymbolicObjectiveFunction',arcarg,repflag);
else
    % return optimal control value evaluated at 'depvar'
    for ii=1:arcn
        arcp=arcpos(1,ii):arcpos(2,ii)+1;
        o=feval(docObj,'ObjectiveFunction',indepvar(arcp(2:end)),depvar(:,arcp),par,arcarg(ii));
        if connectflag
            varargout{1}(arcp(1:end-1))=o;
        else
            varargout{ii}=o;
        end
    end
end
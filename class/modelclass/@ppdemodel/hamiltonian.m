function varargout=hamiltonian(ppdeObj,solObj,varargin)
%

if isempty(ppdeObj)
    varargout{1}=[];
    return
end
if nargin==1
    solObj=[];
end
par=parametervalue(ppdeObj);
if isstruct(solObj)
    try
        arcarg=solObj.arcarg;
        indepvar=solObj.independentvar;
        depvar=solObj.dependentvar;
        arcpos=solObj.arcposition;
        femdat=solObj.femdata;
        arcn=numel(arcarg);
    catch
        ocmaterror(['Field is missing.\n' lasterr])
    end
elseif isppdeprimitive(solObj) || isppdetrajectory(solObj) || isppdeasymptotic(solObj)
    femdat=femdata(solObj);
    arcarg=arcargument(solObj);
    indepvar=time(ppdeObj,solObj,1);
    depvar=dependentvar(solObj);
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
    ocmatmsg('Symbolic representation not implemented yet.')
    return
    %varargout{1}=feval(ppdeObj,'SymbolicCanonicalSystem',arcarg,repflag);
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
        %equationnum=canonicalsystemdimension(ppdeObj,arcarg(ii));
        H=feval(ppdeObj,'Hamiltonian',indepvar(arcp),depvar(:,arcp),par,arcarg(ii),femdat);
        if connectflag
            varargout{1}(:,arcp)=H;
        else
            varargout{ii}=H;
        end
    end
end

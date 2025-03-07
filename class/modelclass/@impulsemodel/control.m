function varargout=control(ocObj,solObj,varargin)
%
% CONTROL returns the control values.
% 
% U=CONTROL(OCOBJ) OCOBJ is a stdocmodel class. U is a cell array of
% strings consisting of the terms derived from the maximum Hamiltonian
% condition (explicit control) or control variable names (implicit
% control).  
% 
% U=CONTROL(OCOBJ,SOLOBJ) SOLOBJ is a solution object
% (dynprimitive,octrajectory,ocasymptotic,occurve) a structure, including 
% at least the fieldnames 'y' and 'arcposition'. If SOLOBJ is empty the
% output is the same as for CONTROL(OCOBJ). Otherwise the control values
% are returned. If SOLOBJ is an octrajectory consisting of multiple arcs U
% is a cell array of matrices, with the control values for each arc
% separately. 
%
% U=CONTROL(OCOBJ,SOLOBJ,CONNECF) CONNECF=0 is the same as
% CONTROL(OCOBJ,SOLOBJ). If CONNECF is 1 multiple arcs are returned
% "connected", i.e., the control values of all arcs are returned in one
% matrix.

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
    indepvar=time(ocObj,solObj,1);
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    arcn=arcnum(solObj);
end

if isempty(solObj)
    % return symbolic expressions
    % if repflag=1 control values/Lagrange multipliers are replaced by the
    % corresponding symbolic expressions for the optimal values
    arcarg=[];
    if nargin>=3
        arcarg=varargin{1};
    end
    if isempty(arcarg)
        arcarg=0;
    end
    varargout{1}=feval(ocObj,'SymbolicOptimalControl',arcarg);
else
    % return value of the canonical system evaluated at 'depvar'
    connectflag=[];
    if nargin>=3
        connectflag=varargin{1};
    end
    if isempty(connectflag)
        connectflag=0;
    end
    ctrlnum=controlnum(ocObj);
    for ii=1:arcn
        arcp=arcpos(1,ii):arcpos(2,ii);
        u=feval(ocObj,'OptimalControl',indepvar(arcp-1),depvar(:,arcp),par,arcarg(ii),[]);
        if connectflag
            varargout{1}(1:ctrlnum,arcp-1)=u;
        else
            varargout{ii}=u;
        end
    end
end

function varargout=hamiltonian(docObj,solObj,varargin)
%
% HAMILTONIAN returns the symbolic Hamiltonian or evaluated at a solution object.
% 
% H=HAMILTONIAN(OCOBJ) OCOBJ is a stdocmodel class. H is a string of the
% formally compounded (extended) Hamiltonian.
% 
% H=HAMILTONIAN(OCOBJ,[],ARCARG,1) ARCARG is the identifier for the used
% combination of active and inactive constraints. Then, H is given by the
% string expression for the maximized Hamiltonian (control variables are
% replaced by its explicit terms). The argument 1 triggers the replacement
% of the control variables (and in case of active constraints of the
% Lagrangian mutlipliers).
% 
% H=HAMILTONIAN(OCOBJ,SOLOBJ) SOLOBJ is a solution object
% (dynprimitive,octrajectory,ocasymptotic,occurve) a structure, including 
% at least the fieldnames 'y' and 'arcposition'. If SOLOBJ is empty the
% output is the same as for HAMILTONIAN(OCOBJ). Otherwise the Hamiltonian
% evaluated at the solution object is returned. If SOLOBJ is an
% octrajectory consisting of multiple arcs H is a cell array of vectors,
% with the Hamiltonian values for each arc separately. 
%
% H=HAMILTONIAN(OCOBJ,SOLOBJ,CONNECF) CONNECF=0 is the same as
% HAMILTONIAN(OCOBJ,SOLOBJ). If CONNECF is 1 multiple arcs are returned
% "connected", i.e., the Hamiltonian values of all arcs are returned in one
% matrix.
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
    varargout{1}=feval(docObj,'SymbolicHamiltonian',arcarg,repflag);
else
    % return optimal control value evaluated at 'depvar'
    for ii=1:arcn
        arcp=arcpos(1,ii):arcpos(2,ii)+1;
        H=feval(docObj,'Hamiltonian',indepvar(arcp(2:end)),depvar(:,arcp),par,arcarg(ii));
        if connectflag
            varargout{1}(arcp(1:end-1))=H;
        else
            varargout{ii}=H;
        end
    end
end
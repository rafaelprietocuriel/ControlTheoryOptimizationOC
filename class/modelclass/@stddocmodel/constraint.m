function varargout=constraint(docObj,solObj,varargin)
%
% CONSTRAINT returns the value of the constraint.
% 
% U=CONSTRAINT(OCOBJ) OCOBJ is a stdocmodel class. U is a cell array of
% strings consisting of the terms derived from the maximum Hamiltonian
% condition (explicit control) or control variable names (implicit
% control).  
% 
% U=CONSTRAINT(OCOBJ,SOLOBJ) SOLOBJ is a solution object
% (dynprimitive,octrajectory,ocasymptotic,occurve) a structure, including 
% at least the fieldnames 'y' and 'arcposition'. If SOLOBJ is empty the
% output is the same as for CONTROL(OCOBJ). Otherwise the control values
% are returned. If SOLOBJ is an octrajectory consisting of multiple arcs U
% is a cell array of matrices, with the control values for each arc
% separately. 
%
% U=CONSTRAINT(OCOBJ,SOLOBJ,CONNECF) CONNECF=0 is the same as
% CONSTRAINT(OCOBJ,SOLOBJ). If CONNECF is 1 multiple arcs are returned
% "connected", i.e., the control values of all arcs are returned in one
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
    arcarg=[];
    if nargin>=3
        arcarg=varargin{1};
    end
    if isempty(arcarg)
        arcarg=0;
    end
    varargout{1}=feval(docObj,'SymbolicConstraint',arcarg);
else
    for ii=1:arcn
        arcp=arcpos(1,ii):arcpos(2,ii)+1;
        c=feval(docObj,'Constraint',indepvar(arcp(2:end)),depvar(:,arcp),par,arcarg(ii));
        if connectflag
            varargout{1}(1:ctrlnum,arcp(1:end-1))=c;
        else
            varargout{ii}=c;
        end
    end
end

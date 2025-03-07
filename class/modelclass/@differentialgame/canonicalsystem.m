function varargout=canonicalsystem(dgObj,solObj,varargin)
%
% CANONICALSYSTEM returns the value of the canonical system
%
% CANONICALSYSTEM(OCOBJ) returns the symbolic expression of the canonical
% system in its general form, i.e., control variables (Lagrangian
% multipliers) are not substituted.
%
% D=CANONICALSYSTEM(OCOBJ,SOLOBJ) SOLOBJ is a solution object
% (dynprimitive,octrajectory,ocasymptotic,occurve) a structure, including 
% at least the fieldnames 'y' and 'arcposition'. If SOLOBJ is empty the
% output is the same as for CANONICALSYSTEM(OCOBJ). Otherwise the values of
% the canonical system are returned. If SOLOBJ is an octrajectory
% consisting of multiple arcs D is a cell array of matrices, with the
% values of the canonical system for each arc separately. 
%
% D=CANONICALSYSTEM(OCOBJ,SOLOBJ,CONNECF) CONNECF=0 is the same as
% CANONICALSYSTEM(OCOBJ,SOLOBJ). If CONNECF is 1 multiple arcs are returned
% "connected", i.e., the values of the canonical system of all arcs are
% returned in one matrix.
%
% D=CANONICALSYSTEM(OCOBJ,[],ARCARG,S) if 'S' is zero the symbolic
% expressions for the control values (Lagrange multipliers) are not
% replaced. Otherwise the values are replaced by its symbolic expressions.
% 

if isempty(dgObj)
    varargout{1}=[];
    return
end
if nargin==1
    solObj=[];
end
par=parametervalue(dgObj);
if isstruct(solObj)
    try
        arcarg=solObj.arcarg;
        indepvar=solObj.independentvar;
        depvar=solObj.dependentvar;
        arcpos=solObj.arcposition;
        arcn=numel(arcarg);
    catch
        ocmaterror(['Field is missing.\n' lasterr])
    end
elseif isdynprimitive(solObj) || isoctrajectory(solObj) || isocasymptotic(solObj)
    arcarg=arcargument(solObj);
    indepvar=time(dgObj,solObj,1);
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    arcn=arcnum(solObj);
elseif  isoccurve(solObj)
    arcarg=arcargument(solObj);
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    if isempty(arcpos)
        arcarg=unique(arcarg);
        arcpos=[1 size(depvar,2)];
        arcn=1;
    else
        arcpos=arcposition(solObj);
        try
            arcn=arcnum(solObj);
        catch
            arcn=1;
        end
    end
    indepvar=zeros(1,size(depvar,2));
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
    varargout{1}=feval(dgObj,'SymbolicCanonicalSystem',arcarg,repflag);
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
        equationnum=canonicalsystemequationnum(dgObj);
        dxdt=feval(dgObj,'CanonicalSystem',indepvar(arcp),depvar(:,arcp),par,arcarg(ii));
        if connectflag
            varargout{1}(1:equationnum,arcp)=dxdt;
        else
            varargout{ii}=dxdt;
        end
    end
end

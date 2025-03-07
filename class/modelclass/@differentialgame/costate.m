function varargout=costate(dgObj,solObj,varargin)
%
% COSTATE returns the costate values/variables.
% 
% L=COSTATE(OCOBJ) OCOBJ is a stdocmodel class. L is a cell array of
% strings consisting of the costate variable names. 
% 
% L=COSTATE(OCOBJ,SOLOBJ) SOLOBJ is a solution object
% (dynprimitive,octrajectory,ocasymptotic,occurve) a structure, including 
% at least the fieldnames 'y' and 'arcposition'. If SOLOBJ is empty the
% output is the same as for COSTATE(OCOBJ). Otherwise the costate values
% are returned. If SOLOBJ is an octrajectory consisting of multiple arcs L
% is a cell array of matrices, with the costate values for each arc
% separately. 
%
% L=COSTATE(OCOBJ,SOLOBJ,CONNECF) CONNECF=0 is the same as
% COSTATE(OCOBJ,SOLOBJ). If CONNECF is 1 multiple arcs are returned
% "connected", i.e., the costate values of all arcs are returned in one
% matrix.

connectflag=[];
if isempty(dgObj)
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
pn=playernum(dgObj);

if isstruct(solObj)
    try
        depvar=solObj.dependentvar;
        arcpos=solObj.arcposition;
    catch
        ocmaterror('If second input argument is a structure the fields ''dependentvar'', ''independentvar'' and ''arcarg'' have to exist!')
    end
elseif isdynprimitive(solObj) || isoctrajectory(solObj) || isocasymptotic(solObj)
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    arcn=arcnum(solObj);
elseif isoccurve(solObj)
    arcn=1;
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    if isempty(arcpos)
        arcpos=[1 size(depvar,2)].';
    end
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
    depvar=dependentvariable(dgObj,arcarg);
    varargout{1}=depvar(costatecoord(dgObj));
else
    for ii=1:arcn
        arcp=arcpos(1,ii):arcpos(2,ii);
        lcoord=costatecoordinate(dgObj);
        l=depvar(lcoord,arcp);
        if connectflag
            varargout{1}(:,arcp)=l;
        else
            varargout{ii}=l;
        end
    end
end
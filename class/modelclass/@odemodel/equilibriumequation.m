function varargout=equilibriumequation(odeObj,solObj,varargin)
%
% EQUILIBRIUMEQUATION returns the value of the canonical system
%
% F=EQUILIBRIUMEQUATION(odeObj,SOLOBJ) SOLOBJ is a solution object
% (dynprimitive, octrajectory, ocasymptotic, occurve) a structure,
% including  the mandatory fieldnames 'y' and 'arcposition'. The values of
% the equilibrium equations are returned. 

if isempty(odeObj)
    varargout{1}=[];
    return
end
if nargin==1
    solObj=[];
end
if isstruct(solObj)
    try
        arcarg=solObj.arcarg;
        depvar=solObj.dependentvar;
        arcpos=solObj.arcposition;
        arcn=numel(arcarg);
    catch
        ocmaterror(['Field is missing.\n' lasterr])
    end
elseif isdynprimitive(solObj) || isoctrajectory(solObj) || isocasymptotic(solObj)
    arcarg=arcargument(solObj);
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
        arcn=arcnum(solObj);
    end
end

if ~isempty(solObj)
    % returns the value of the e evaluated at 'depvar'
    for ii=1:arcn
        arcp=arcpos(1,ii):arcpos(2,ii);
        dxdt=feval(odeObj,'EquilibriumEquation',depvar(:,arcp),parametervalue(odeObj),arcarg(ii));
        varargout{ii}=dxdt;
    end
else
    arcarg=[];
    if nargin>=3
        arcarg=varargin{1};
    end
    if isempty(arcarg)
        arcarg=0;
    end
    varargout{1}=feval(odeObj,'SymbolicEquilibriumEquation',arcarg);
end

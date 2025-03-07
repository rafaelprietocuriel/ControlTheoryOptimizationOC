function varargout=equilibriumequation(ppdeObj,solObj,varargin)
%
% EQUILIBRIUMEQUATION returns the value of the canonical system
%
% F=EQUILIBRIUMEQUATION(OCOBJ,SOLOBJ) SOLOBJ is a solution object
% (dynprimitive, octrajectory, ocasymptotic, occurve) a structure,
% including  the mandatory fieldnames 'y' and 'arcposition'. The values of
% the equilibrium equations are returned. 

if isempty(ppdeObj)
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
        femdat=solObj.femdata;
        arcn=numel(arcarg);
    catch
        ocmaterror(['Field is missing.\n' lasterr])
    end
elseif ispdeprimitive(solObj) || ispdetrajectory(solObj) || ispdeasymptotic(solObj)
    arcarg=arcargument(solObj);
    depvar=dependentvar(solObj);
    arcpos=arcposition(solObj);
    femdat=femdata(solObj);
    arcn=arcnum(solObj);
end

% returns the value of the e evaluated at 'depvar'
for ii=1:arcn
    arcp=arcpos(1,ii):arcpos(2,ii);
    dxdt=feval(ppdeObj,'EquilibriumEquation',depvar(:,arcp),parametervalue(ppdeObj),arcarg(ii),femdat);
    varargout{ii}=dxdt;
end

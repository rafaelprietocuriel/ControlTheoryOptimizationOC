function varargout=equilibriumequation(ocObj,solObj,varargin)
%
% EQUILIBRIUMEQUATION returns the value of the canonical system
%
% F=EQUILIBRIUMEQUATION(OCOBJ,SOLOBJ) SOLOBJ is a solution object
% (dynprimitive, octrajectory, ocasymptotic, occurve) a structure,
% including  the mandatory fieldnames 'y' and 'arcposition'. The values of
% the equilibrium equations are returned. 
equationfile='';
symbolicequationfile='';

if isempty(ocObj)
    varargout{1}=[];
    return
end
if nargin==1
    solObj=[];
end
equationfileidx=find(strcmpi(varargin,'equationfile'));
symbolicequationfileidx=find(strcmpi(varargin,'symbolicequationfile'));
if ~isempty(symbolicequationfileidx)
    symbolicequationfile=varargin{symbolicequationfileidx+1};
end
if ~isempty(equationfileidx)
    equationfile=varargin{equationfileidx+1};
end
if isempty(equationfile)
    equationfile='CanonicalSystem';
end
if isempty(symbolicequationfile)
    symbolicequationfile='SymbolicCanonicalSystem';
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
        if strcmp(equationfile,'CanonicalSystem')
            dxdt=feval(ocObj,equationfile,0,depvar(:,arcp),parametervalue(ocObj),arcarg(ii));
        else
            dxdt=feval(ocObj,equationfile,depvar(:,arcp),parametervalue(ocObj),arcarg(ii));
        end
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
    if strcmp(symbolicequationfile,'SymbolicCanonicalSystem')
        varargout{1}=feval(ocObj,symbolicequationfile,arcarg,1);
    else
        varargout{1}=feval(ocObj,symbolicequationfile,arcarg);
    end
end

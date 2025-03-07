function sol=initstatgrad_ES(ocObj,ocgPt,contpar,varargin)
clear global OCGRADSOL OCGRADCONT
global OCGRADSOL OCGRADCONT
targetvalue=[];
freeparametervalue=[];

targetvalueidx=find(strcmpi(varargin,'targetvalue'));
freeparametervalueidx=find(strcmpi(varargin,'freeparametervalue'));
if ~isempty(freeparametervalueidx)
    freeparametervalue=varargin{freeparametervalueidx+1};
end
if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end
if nargin==3
    targetvalue=[];
end
if isnumeric(ocgPt)
    ocgPt0.variable.y=ocgPt;
    ocgPt0.initargument=[];
    ocgPt0.argument=[];
    ocgPt0.modelparameter=parametervalue(ocObj);
    ocgPt0.modelname=modelname(ocObj);
    ocgPt=ocgradtrajectory(ocgPt0);
end

modelfunction=modelspecificfunc(ocObj,'4GradContinuation');
par=parametervalue(ocObj);
func=modelfunction();

OCGRADSOL.objective=func{1}{1};
OCGRADSOL.inequalityconstraint=func{1}{2};
OCGRADSOL.equalityconstraint=func{1}{3};
OCGRADSOL.gradientobjective=func{2}{1};
OCGRADSOL.gradientic=func{2}{2};
OCGRADSOL.gradientec=func{2}{3};
OCGRADSOL.qpsolution=func{3};
OCGRADSOL.nesterovsolution=func{4};

OCGRADSOL.plotcontinuation=func{18};
OCGRADSOL.datapath=func{20};
OCGRADSOL.saveintermediatefiles=func{21};

OCGRADSOL.statenum=statenum(ocObj);
OCGRADSOL.equalitynum=equalityconstraintnum(ocObj);
OCGRADSOL.inequalitynum=inequalityconstraintnum(ocObj);
if ischar(contpar)
    contpar=parameterindex(ocObj,contpar);
end
if isempty(freeparametervalue)
    freeparametervalue=par(contpar);
else
    par(contpar)=freeparametervalue;
end
OCGRADSOL.modelparameter=par;
OCGRADSOL.modelname=modelname(ocObj);


sol.extremal=generatestaticextremal(ocgPt);
sol.parameters=freeparametervalue;
OCGRADSOL.continuationindex=1;

if isnumeric(contpar)
    OCGRADSOL.contparidx=contpar;
end

%OCGRADSOL.equalityindex=1:p;
%OCGRADSOL.d=ones(n,1);

OCGRADCONT.targetvalue=targetvalue;

OCGRADSOL.d=eye(OCGRADSOL.statenum);
pathname=OCGRADSOL.datapath();
[resultfile,globalvarfile]=OCGRADSOL.saveintermediatefiles();
OCGRADSOL.basicresultfilename=fullocmatfile(pathname,resultfile);
OCGRADSOL.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCGRADCONT.basicmoviefilename=fullocmatfile(pathname,'AnimCont');


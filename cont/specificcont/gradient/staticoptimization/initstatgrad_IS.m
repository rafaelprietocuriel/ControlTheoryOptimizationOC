function sol=initstatgrad_IS(ocObj,ocgMTrj,contpar,varargin)
clear global OCGRADCONT OCGRADSOL
global OCGRADCONT OCGRADSOL
freeparameter=[];
targetvalue=[];
targetindex=[];

freeparameteridx=find(strcmpi(varargin,'freeparameter')); %
targetvalueidx=find(strcmpi(varargin,'targetvalue'));
targetindexidx=find(strcmpi(varargin,'targetindex'));
if ~isempty(freeparameteridx)
    freeparameter=varargin{freeparameteridx+1};
end
if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end
if ~isempty(targetindexidx)
    targetindex=varargin{targetindexidx+1};
end
if ischar(targetindex)
    targetindex=parameterindex(ocObj,targetindex);
end
if isempty(freeparameter)
    sol=[];
    return
end
ocgMTrj=ocgradmultipath(ocgMTrj);

if ~strcmp(modelname(ocObj),modelname(ocgMTrj))
    ocmatmsg('Warning, model name of extremal and model are different')
end

OCGRADCONT.modelfunc=modelspecificfunc(ocObj,'4GradContinuation');
funch=OCGRADCONT.modelfunc(); % model specific function handles for saddle path continuation

OCGRADSOL.objective=funch{1}{1};
OCGRADSOL.inequalityconstraint=funch{1}{2};
OCGRADSOL.equalityconstraint=funch{1}{3};
OCGRADSOL.gradientobjective=funch{2}{1};
OCGRADSOL.gradientic=funch{2}{2};
OCGRADSOL.gradientec=funch{2}{3};
OCGRADSOL.gradientparobjective=funch{2}{4};
OCGRADSOL.qpsolution=funch{3};
OCGRADSOL.nesterovsolution=funch{4};

OCGRADSOL.plotcontinuation=funch{18};
OCGRADSOL.datapath=funch{20};
OCGRADSOL.saveintermediatefiles=funch{21};

% filenames for storing intermediate continuation results
pathname=OCGRADSOL.datapath();
[resultfile,globalvarfile]=OCGRADSOL.saveintermediatefiles();
OCGRADSOL.basicresultfilename=fullocmatfile(pathname,resultfile);
OCGRADSOL.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);

OCGRADCONT.modelname=modelname(ocObj);
OCGRADCONT.modelparameter=parametervalue(ocObj);

OCGRADSOL.degree=degree(ocgMTrj);
OCGRADCONT.octype=octype(ocgMTrj);
OCGRADCONT.targetvalue=targetvalue;

contindex=parameterindex(ocObj,contpar);
OCGRADCONT.contindex=contindex;

OCGRADCONT.freeparameterindex=parameterindex(ocObj,freeparameter);

sol.extremal=generategradextremal(ocgMTrj);
sol.parameters=parametervalue(ocObj,freeparameter);
OCGRADCONT.freeparametercoordinate=1:length(sol.parameters);

sol.parameters=[sol.parameters;parametervalue(ocObj,contindex)];
OCGRADCONT.continuationcoordinate=length(sol.parameters);

OCGRADSOL.statenum=statenum(ocObj);
OCGRADSOL.equalitynum=equalityconstraintnum(ocObj);
OCGRADSOL.inequalitynum=inequalityconstraintnum(ocObj);

%OCGRADCONT.equalityindex=1:p;
OCGRADCONT.d=eye(OCGRADSOL.statenum);

if ~isempty(targetvalue) && isempty(targetindex)
    targetindex=contindex;
end
OCGRADCONT.targetindex=targetindex;

OCGRADCONT.statecoordinate{1}=(1:OCGRADSOL.statenum).';
for ii=2:OCGRADSOL.degree
    OCGRADCONT.statecoordinate{ii}=OCGRADCONT.statecoordinate{ii-1}(end)+(1:OCGRADSOL.statenum).';
end



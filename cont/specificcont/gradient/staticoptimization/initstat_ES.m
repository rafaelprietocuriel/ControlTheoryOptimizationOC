function sol=initstat_ES(ocObj,ocPt,contpar,varargin)
clear global OCSTATSOL OCSTATCONT
global OCSTATSOL OCSTATCONT
targetvalue=[];
freeparametervalue=[];
targettype='';
targetcoordinate=[];

targetvalueidx=find(strcmpi(varargin,'targetvalue'));
targettypeidx=find(strcmpi(varargin,'targettype'));
targetcoordinateidx=find(strcmpi(varargin,'targetcoordinate'));
freeparametervalueidx=find(strcmpi(varargin,'freeparametervalue'));
if ~isempty(freeparametervalueidx)
    freeparametervalue=varargin{freeparametervalueidx+1};
end
if ~isempty(targetvalueidx)
    targetvalue=varargin{targetvalueidx+1};
end
if ~isempty(targettypeidx)
    targettype=varargin{targettypeidx+1};
end
if ~isempty(targetcoordinateidx)
    targetcoordinate=varargin{targetcoordinateidx+1};
end

if ischar(contpar)
    contpar=parameterindex(ocObj,contpar);
end
if isempty(targettype) && ~isempty(targetvalue)
    targettype='parameter';
    targetcoordinate=contpar;
end
n=statenum(ocObj);
m=constraintnum(ocObj);
if isnumeric(ocPt)
    if length(ocPt)==n
        ocPt0.x=ocPt;
        ocPt0.lm=zeros(m,1);
    elseif length(ocPt)==n+m
        ocPt0.x=ocPt(1:n);
        ocPt0.lm=ocPt(n+1:m+n);
    else
        sol=[];
        return
    end
    ocPt0.modelparameter=parametervalue(ocObj);
    ocPt0.modelname=modelname(ocObj);
    ocPt=staticextremal(ocPt0);
end

modelfunction=modelspecificfunc(ocObj,'4Continuation');
par=parametervalue(ocObj);
func=modelfunction();

OCSTATCONT.targettype=targettype;
OCSTATCONT.targetcoordinate=targetcoordinate;

OCSTATSOL.firstordernc=func{1};
OCSTATSOL.jacobianfirstordernc=func{2};
OCSTATSOL.parameterjacobianfirstordernc=func{3};
OCSTATSOL.lagrangefunction=func{5};
OCSTATSOL.targetfunction=func{19};
OCSTATSOL.plotcontinuation=func{18};
OCSTATSOL.datapath=func{20};
OCSTATSOL.saveintermediatefiles=func{21};

OCSTATSOL.statenum=n;
OCSTATSOL.constraintnum=m;
if isempty(freeparametervalue)
    freeparametervalue=par(contpar);
else
    par(contpar)=freeparametervalue;
end
OCSTATSOL.modelparameter=par;
OCSTATSOL.modelname=modelname(ocObj);

OCSTATSOL.variableindex=1:n;
OCSTATSOL.lagrangemulitplierindex=n+1:n+m;

sol.coeff=[state(ocPt);lagrangemultiplier(ocPt);freeparametervalue];
OCSTATSOL.continuationindex=n+m+1;

if isnumeric(contpar)
    OCSTATSOL.parametercoordinate=contpar;
end

%OCSTATSOL.equalityindex=1:p;
%OCSTATSOL.d=ones(n,1);

OCSTATCONT.targetvalue=targetvalue;

pathname=OCSTATSOL.datapath();
[resultfile,globalvarfile]=OCSTATSOL.saveintermediatefiles();
OCSTATSOL.basicresultfilename=fullocmatfile(pathname,resultfile);
OCSTATSOL.basicglobalvarfilename=fullocmatfile(pathname,globalvarfile);
OCSTATCONT.basicmoviefilename=fullocmatfile(pathname,'AnimCont');


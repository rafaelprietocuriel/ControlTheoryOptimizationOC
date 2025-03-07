function varargout=lagrangemultiplier(ocgTrj,varargin)
%
%
if nargin==1
    ocObj=stdocmodel(modelname(ocgTrj),[],[],0);
end

if nargin>=2
    ocObj=varargin{1};
end
par=modelparameter(ocgTrj);
ocObj=changeparametervalue(ocObj,par);

arcarg=arcargument(ocgTrj);
indepvar=time(ocObj,ocgTrj,1);
depvar=dependentvar(ocgTrj);
arcpos=arcposition(ocgTrj);
arcn=arcnum(ocgTrj);

out=cell(2,arcn);
for ii=1:arcn
    [lmmc,lmsc]=feval(ocObj,'LagrangeMultiplier',indepvar(arcpos(1,ii):arcpos(2,ii)),depvar{ii},par,arcarg(ii));
    out{1}{ii}=lmmc;
    out{2}{ii}=lmsc;
end
varargout{1}=out{1};
varargout{2}=out{2};

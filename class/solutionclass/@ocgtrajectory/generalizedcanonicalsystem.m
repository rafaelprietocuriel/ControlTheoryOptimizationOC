function out=generalizedcanonicalsystem(ocgTrj,varargin)
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

n=statenum(ocObj);
m=controlnum(ocObj);
out=cell(1,arcn);
for ii=1:arcn
    coord=implicitcontrolcoordinate(ocObj,arcarg(ii));
    dxdt=feval(ocObj,'CanonicalSystemGeneralized',indepvar(arcpos(1,ii):arcpos(2,ii)),depvar{ii},par,arcarg(ii));
    out{ii}=zeros(2*n+m,size(dxdt,2));
    out{ii}([1:2*n 2*n+coord],:)=dxdt;
end

function out=equilibriumequation(ocgTrj,varargin)
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
depvar=dependentvar(ocgTrj);
arcn=arcnum(ocgTrj);

out=cell(1,arcn);
for ii=1:arcn
    out{ii}=feval(ocObj,'EquilibriumEquation',depvar{ii},par,arcarg(ii));
end

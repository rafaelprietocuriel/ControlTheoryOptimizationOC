function out=dhamiltoniandu(ocgTrj,varargin)
%
%
ocObj=loadmodel(ocgTrj);
if isempty(ocObj)
    out=[];
    return
end
par=modelparameter(ocgTrj);
ocObj=changeparametervalue(ocObj,par);

arcarg=arcargument(ocgTrj);
indepvar=time(ocObj,ocgTrj,1);
depvar=dependentvar(ocgTrj);
arcpos=arcposition(ocgTrj);
arcn=arcnum(ocgTrj);

out=cell(1,arcn);
for ii=1:arcn
    out{ii}=feval(ocObj,'DHamiltonianDU',indepvar(arcpos(1,ii):arcpos(2,ii)),depvar{ii},par,arcarg(ii));
end

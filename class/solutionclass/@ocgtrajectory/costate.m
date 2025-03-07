function out=costate(ocgTrj,varargin)
%
%
ocObj=loadmodel(ocgTrj);
if isempty(ocObj)
    out=[];
    return
end
par=modelparameter(ocgTrj);
ocObj=changeparametervalue(ocObj,par);
n=statenum(ocObj);

depvar=dependentvar(ocgTrj);
arcn=arcnum(ocgTrj);

out=cell(1,arcn);
for ii=1:arcn
    out{ii}=depvar{ii}(n+1:2*n,:);
end

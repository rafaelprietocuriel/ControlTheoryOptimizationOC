function out=derivativegeneralizedcontrol(ocgTrj)
%
% DERIVATIVEGENERALIZEDCONTROL returns the derivative of the generalized
% control(s) (control(s) and Lagrangemultiplier with respect to the
% state(s) and costate(s).
%
ocObj=loadmodel(ocgTrj);
if isempty(ocObj)
    ocmatmsg('Input model is empty.\n')
    out=[];
    return
end

if isempty(ocgTrj)
    ocmatmsg('Input solution is empty.\n')
    out=[];
    return
end
par=modelparameter(ocgTrj);
ocObj=changeparametervalue(ocObj,par);

arcarg=arcargument(ocgTrj);
t=time(ocgTrj);
depvar=dependentvar(ocgTrj);
arcpos=arcposition(ocgTrj);
arcn=arcnum(ocgTrj);

out=cell(1,arcn);
for ii=1:arcn
    out{ii}=feval(ocObj,'DImplicitControlDX',t(arcpos(1,ii):arcpos(2,ii)),depvar{ii},par,arcarg);
end


function out=daeorder(ocTrj)
out=[];

if isempty(ocTrj)
    return
end

if isfield(ocTrj.solverinfo,'daeorder')
    out=ocTrj.solverinfo.daeorder;
else
    ocmatmsg('Warning: order is calculated from model properties.\n')
    ocObj=stdocmodel(modelname(ocTrj));
    out=[ones(2*statenum(ocObj),1);zeros(controlnum(ocObj)+inequalitycontrolconstraintnum(ocObj),1)];
end
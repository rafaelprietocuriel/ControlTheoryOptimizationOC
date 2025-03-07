function out=inftimetransformation(ocTrj)
out=[];
if isfield(ocTrj.solverinfo,'inftimetransformation')
    out=ocTrj.solverinfo.inftimetransformation;
end
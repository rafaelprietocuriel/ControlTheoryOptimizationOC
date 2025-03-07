function out=pathtype(ocTrj)

out=[];
if isfield(ocTrj.solverinfo,'pathtype')
    out=ocTrj.solverinfo.pathtype;
end
function out=interpolationmethod(ocTrj)
out=[];

if isempty(ocTrj)
    return
end

if isfield(ocTrj.solverinfo,'interpolationmethod')
    out=ocTrj.solverinfo.interpolationmethod;
end
function out=collocationmethod(ocTrj)
out=[];

if isempty(ocTrj)
    return
end

if isfield(ocTrj.solverinfo,'collocationmethod')
    out=ocTrj.solverinfo.collocationmethod;
end
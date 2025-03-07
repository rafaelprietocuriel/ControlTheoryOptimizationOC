function out=method(ocTrj)
out=[];

if isempty(ocTrj)
    return
end

if isfield(ocTrj.solverinfo,'method')
    out=ocTrj.solverinfo.method;
end
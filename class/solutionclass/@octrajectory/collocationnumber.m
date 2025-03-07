function out=collocationnumber(ocTrj)
out=[];

if isempty(ocTrj)
    return
end

if isfield(ocTrj.solverinfo,'collocationnumber')
    out=ocTrj.solverinfo.collocationnumber;
end
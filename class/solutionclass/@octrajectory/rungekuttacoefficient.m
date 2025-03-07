function out=rungekuttacoefficient(ocTrj)
out=[];

if isempty(ocTrj)
    return
end

if isfield(ocTrj.solverinfo,'coeff')
    out=ocTrj.solverinfo.coeff;
end
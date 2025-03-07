function out=conttype(ocTrj)


if ~isfield(ocTrj.solverinfo,'conttype')
    out='unknown';
    return
end
out=ocTrj.solverinfo.conttype;
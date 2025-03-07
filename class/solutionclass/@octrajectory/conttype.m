function out=conttype(ocTrj)

if ~isfield(ocTrj.solverinfo,'conttype') && ~isfield(ocTrj.solverinfo,'contclass')
    out='unknown';
    return
end
if isfield(ocTrj.solverinfo,'conttype')
    out=ocTrj.solverinfo.conttype;
else
    out=ocTrj.solverinfo.contclass;
end
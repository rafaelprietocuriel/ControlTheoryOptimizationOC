function out=contclass(ocTrj)

if ~isfield(ocTrj.solverinfo,'contclass')
    out='unknown';
    return
end
out=ocTrj.solverinfo.contclass;

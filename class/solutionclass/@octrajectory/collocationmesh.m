function out=collocationmesh(ocTrj)
out=[];

if isempty(ocTrj)
    return
end

if isfield(ocTrj.solverinfo,'tcolmesh')
    out=ocTrj.solverinfo.tcolmesh;
end
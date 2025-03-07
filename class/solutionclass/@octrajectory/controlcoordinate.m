function out=controlcoordinate(ocTrj)
out=[];

if isempty(ocTrj)
    return
end

if isfield(ocTrj.solverinfo,'controlcoordinate')
    out=ocTrj.solverinfo.controlcoordinate;
end
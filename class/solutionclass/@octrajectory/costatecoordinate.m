function out=costatecoordinate(ocTrj)
out=[];

if isempty(ocTrj)
    return
end

if isfield(ocTrj.solverinfo,'costatecoordinate')
    out=ocTrj.solverinfo.costatecoordinate;
end
function out=statecoordinate(ocTrj)
out=[];

if isempty(ocTrj)
    return
end

if isfield(ocTrj.solverinfo,'statecoordinate')
    out=ocTrj.solverinfo.statecoordinate;
end
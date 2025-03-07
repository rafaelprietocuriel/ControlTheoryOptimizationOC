function out=objectivevaluecoordinate(ocTrj)
out=[];

if isempty(ocTrj)
    return
end

if isfield(ocTrj.solverinfo,'objectivevaluecoordinate')
    out=ocTrj.solverinfo.objectivevaluecoordinate;
elseif isfield(ocTrj.solverinfo,'objectivevaluecoord')
    out=ocTrj.solverinfo.objectivevaluecoord;
end
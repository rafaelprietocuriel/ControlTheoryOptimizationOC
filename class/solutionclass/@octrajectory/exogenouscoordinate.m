function out=exogenouscoordinate(ocTrj)
out=[];

if isempty(ocTrj)
    return
end

if isfield(ocTrj.solverinfo,'exogenousdynamicscoordinate')
    out=ocTrj.solverinfo.exogenousdynamicscoordinate;
end
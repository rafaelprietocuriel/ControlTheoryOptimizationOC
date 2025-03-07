function out=lagrangemultipliercoordinate(ocTrj)
out=[];

if isempty(ocTrj)
    return
end

if isfield(ocTrj.solverinfo,'lagrangemultipliercoordinate')
    out=ocTrj.solverinfo.lagrangemultipliercoordinate;
end
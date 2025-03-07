function out=reducedabnormallagrangemultiplier(ocObj,arcid)
%
if isempty(ocObj)
    out=[];
    return
end
if nargin<2
    arcid=[];
end
if isempty(arcid)
    arcid=0;
end
out=feval(ocObj,'ReducedAbnormalSymbolicLagrangeMultiplier',arcid);

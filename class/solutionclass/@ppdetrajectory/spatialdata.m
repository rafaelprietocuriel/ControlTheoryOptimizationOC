function [spatialx spatialy]=spatialdata(solObj,varargin)
%
spatialx=[];
spatialy=[];
plotflag=[];

if isempty(solObj)
    return
end

if nargin>=2
    plotflag=varargin{1};
end
if isempty(plotflag)
    plotflag=0;
end
if plotflag
    plotflag=1;
end
if isppdeprimitive(solObj) ||  isppdetrajectory(solObj) || isppdeasymptotic(solObj)
    pt=points(solObj);
    tr=triangles(solObj);
end
it1=tr(1,:);
it2=tr(2,:);
it3=tr(3,:);

if plotflag
    spatialx=[pt(1,it1).' pt(1,it2).' pt(1,it3).' pt(1,it1).' NaN*ones(size(it1.'))].';
    spatialy=[pt(2,it1).' pt(2,it2).' pt(2,it3).' pt(2,it1).' NaN*ones(size(it1.'))].';
else
    spatialx=[pt(1,it1).' pt(1,it2).' pt(1,it3).' pt(1,it1).'].';
    spatialy=[pt(2,it1).' pt(2,it2).' pt(2,it3).' pt(2,it1).'].';
end
spatialx=spatialx(:);
spatialy=spatialy(:);
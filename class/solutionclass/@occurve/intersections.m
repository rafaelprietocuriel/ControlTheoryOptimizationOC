function [x0,y0,iout,jout]=intersections(ocCurve1,ocCurve2,coord,varargin)

if nargin==2
    coord=1:2;
end

[x0,y0,iout,jout]=intersections(ocCurve1.y(coord(1),:),ocCurve1.y(coord(2),:),ocCurve2.y(coord(1),:),ocCurve2.y(coord(2),:),varargin{:});
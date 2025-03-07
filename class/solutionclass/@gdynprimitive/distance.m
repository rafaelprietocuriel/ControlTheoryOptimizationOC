function out=distance(gdynPrim1,gdynPrim2,varargin)
x1=state(gdynPrim1);
x2=state(gdynPrim2);
l1=costate(gdynPrim1);
l2=costate(gdynPrim2);
ctrl1=control(gdynPrim1);
ctrl2=control(gdynPrim2);
out=norm([x1{1};l1{1};ctrl1{1}]-[x2{1};l2{1};ctrl2{1}],varargin{:});
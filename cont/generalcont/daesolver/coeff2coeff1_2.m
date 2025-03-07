function [tcolmesh1_2,coeff1_2]=coeff2coeff1_2(tcolmesh,coeff)
global OCMATCONT
% returns the coefficients for MESHDATA1_2 correspoing to coeff for
% MESHDATA

sol.y=coeff2points(tcolmesh,coeff,'grid');
sol.x=tcolmesh(OCMATCONT.MESHDATA.tmeshidx);
tcolmesh1_2=makefinemesh(tcolmesh);
t1_2=tcolmesh1_2(OCMATCONT.MESHDATA1_2.tmeshidx);
coeff1_2=points2coeff(sol,t1_2);
coeff1_2(OCMATCONT.MESHDATA1_2.continuationindex)=coeff(OCMATCONT.MESHDATA.continuationindex);
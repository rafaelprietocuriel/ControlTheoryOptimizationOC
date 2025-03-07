function vec2 = calcmultilinear2(funch,q1,q2,tmesh,coeff,tangent,increment)
if (q1==q2)

    vec2 = Bvv(funch,q1,tmesh,coeff,tangent,increment);

else

    part1 = Bvv(funch,q1+q2,tmesh,coeff,tangent,increment);

    part2 = Bvv(funch,q1-q2,tmesh,coeff,tangent,increment);

    vec2 = (part1-part2)/4.0;

end

function tempvec = Bvv(funch,vq,tmesh,coeff,tangent,increment)
global OCMATCONT


vs=vq/norm(vq);
coeff1=coeff;
coeff1(1:OCMATCONT.HE.numdvariablesmc)=coeff1(1:OCMATCONT.HE.numdvariablesmc)+increment*vs;
coeff2=coeff;
coeff2(1:OCMATCONT.HE.numdvariablesmc)=coeff2(1:OCMATCONT.HE.numdvariablesmc)-increment*vs;
f0 = funch(tmesh,coeff,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun);
f1 = funch(tmesh,coeff1,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun);
f2 = funch(tmesh,coeff2,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun);

tempvec = (f1+f2-2.0*f0)*(norm(vq)/increment)^2;
tempvec=tempvec(1:OCMATCONT.HE.numdvariablesmc);


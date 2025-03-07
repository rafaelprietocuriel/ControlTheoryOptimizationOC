function vec3 = calcmultilinear3(funch,q1,q2,q3,tmesh,coeff,tangent,increment)
if (q1==q2)

    if (q1==q3)

        vec3 = Cvvv(funch,q1,tmesh,coeff,tangent,increment);

    else

        part1 = Cvvv(funch,q1+q3,tmesh,coeff,tangent,increment);

        part2 = Cvvv(funch,q1-q3,tmesh,coeff,tangent,increment);

        part3 = Cvvv(funch,q3,tmesh,coeff,tangent,increment);

        vec3 = (part1 - part2 - 2.0*part3)/6.0;

    end

else

    part1 = Cvvv(funch,q1+q2+q3,tmesh,coeff,tangent,increment);

    part2 = Cvvv(funch,q1+q2-q3,tmesh,coeff,tangent,increment);

    part3 = Cvvv(funch,q1-q2+q3,tmesh,coeff,tangent,increment);

    part4 = Cvvv(funch,q1-q2-q3,tmesh,coeff,tangent,increment);

    vec3 = (part1 - part2 - part3 + part4)/24.0;

end


%----------------------------------------------------

function tempvec = Cvvv(funch,vq,tmesh,coeff,tangent,increment)
global OCMATCONT

vs=vq/norm(vq);
coeff1=coeff;
coeff1(1:OCMATCONT.HE.numdvariablesmc)=coeff1(1:OCMATCONT.HE.numdvariablesmc)+3.0*increment*vs;
coeff2=coeff;
coeff2(1:OCMATCONT.HE.numdvariablesmc)=coeff2(1:OCMATCONT.HE.numdvariablesmc)+increment*vs;
coeff3=coeff;
coeff3(1:OCMATCONT.HE.numdvariablesmc)=coeff3(1:OCMATCONT.HE.numdvariablesmc)-increment*vs;
coeff4=coeff;
coeff4(1:OCMATCONT.HE.numdvariablesmc)=coeff4(1:OCMATCONT.HE.numdvariablesmc)-3.0*increment*vs;


f1 = funch(tmesh,coeff1,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun);
f2 = funch(tmesh,coeff2,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun);
f3 = funch(tmesh,coeff3,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun);
f4 = funch(tmesh,coeff4,tangent,OCMATCONT.ode,OCMATCONT.bc,OCMATCONT.icfun);

tempvec = (f1 - 3.0*f2 + 3.0*f3 - f4)/8.0*(norm(vq)/increment)^3;
tempvec=tempvec(1:OCMATCONT.HE.numdvariablesmc);


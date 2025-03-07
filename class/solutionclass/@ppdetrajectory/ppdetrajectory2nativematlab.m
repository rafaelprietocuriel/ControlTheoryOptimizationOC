function sol=ppdetrajectory2nativematlab(ppdeTrj)
%
% OCTRAJECTORY2NATIVEMATLAB transformation to an ode solution structure
%
% SOL=OCTRAJECTORY2NATIVEMATLAB(OCTRJ) the instance OCTRJ of an
% octrajectory is transformed into a usual MATLAB ode solution structure
% SOL. This specifically concerns the % data used for interpolation, stored
% in the 'discretizationinfo' field of OCTRJ.

sol.x=ppdeTrj.t;
sol.y=ppdeTrj.y;
if isfield(ppdeTrj.discretizationinfo,'parameters')
    sol.parameters=ppdeTrj.discretizationinfo.parameters;
end
switch ppdeTrj.discretizationinfo.method
    case 'bvp4c'
        if isfield(ppdeTrj.discretizationinfo,'yp')
            sol.solver=ppdeTrj.discretizationinfo.method;
            sol.yp=ppdeTrj.discretizationinfo.yp;
        end
    case 'bvp5c'
        if isfield(ppdeTrj.discretizationinfo,'yp') &&  isfield(ppdeTrj.discretizationinfo,'ypmid')
            sol.solver=ppdeTrj.discretizationinfo.method;
            sol.idata.yp=ppdeTrj.discretizationinfo.yp;
            sol.idata.ypmid=ppdeTrj.discretizationinfo.ypmid;
        end
    case 'bvp6c'
        if isfield(ppdeTrj.discretizationinfo,'yp') &&  isfield(ppdeTrj.discretizationinfo,'ypmid')
            sol.solver=ppdeTrj.discretizationinfo.method;
            sol.yp=ppdeTrj.discretizationinfo.yp;
            sol.ypmid=ppdeTrj.discretizationinfo.ypmid;
        end
end
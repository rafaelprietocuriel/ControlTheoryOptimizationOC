function sol=ocmultipath2nativematlab(ocMultiPath)
%
% OCTRAJECTORY2NATIVEMATLAB transformation to an ode solution structure
%
% SOL=OCTRAJECTORY2NATIVEMATLAB(OCTRJ) the instance OCTRJ of an
% octrajectory is transformed into a usual MATLAB ode solution structure
% SOL. This specifically concerns the % data used for interpolation, stored
% in the 'solverinfo' field of OCTRJ.

sol.x=ocMultiPath.x;
sol.y=ocMultiPath.y;
if isfield(ocMultiPath.solverinfo,'parameters')
    sol.parameters=ocMultiPath.solverinfo.parameters;
end
switch ocMultiPath.solver
    case 'bvp4c'
        if isfield(ocMultiPath.solverinfo,'yp')
            sol.solver=ocMultiPath.solver;
            sol.yp=ocMultiPath.solverinfo.yp;
        end
    case 'bvp5c'
        if isfield(ocMultiPath.solverinfo,'yp') &&  isfield(ocMultiPath.solverinfo,'ypmid')
            sol.x=ocMultiPath.x;
            sol.y=ocMultiPath.y;
            sol.solver=ocMultiPath.solver;
            sol.idata.yp=ocMultiPath.solverinfo.yp;
            sol.idata.ypmid=ocMultiPath.solverinfo.ypmid;
        end
    case 'bvp6c'
        if isfield(ocMultiPath.solverinfo,'yp') &&  isfield(ocMultiPath.solverinfo,'ypmid')
            sol.x=ocMultiPath.x;
            sol.y=ocMultiPath.y;
            sol.solver=ocMultiPath.solver;
            sol.yp=ocMultiPath.solverinfo.yp;
            sol.ypmid=ocMultiPath.solverinfo.ypmid;
        end
end
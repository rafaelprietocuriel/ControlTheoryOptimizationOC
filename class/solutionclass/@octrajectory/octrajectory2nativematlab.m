function sol=octrajectory2nativematlab(ocTrj)
%
% OCTRAJECTORY2NATIVEMATLAB transformation to an ode solution structure
%
% SOL=OCTRAJECTORY2NATIVEMATLAB(OCTRJ) the instance OCTRJ of an
% octrajectory is transformed into a usual MATLAB ode solution structure
% SOL. This specifically concerns the % data used for interpolation, stored
% in the 'solverinfo' field of OCTRJ.

sol.x=ocTrj.x;
sol.y=ocTrj.y;
if isfield(ocTrj.solverinfo,'parameters')
    sol.parameters=ocTrj.solverinfo.parameters;
end
switch ocTrj.solver
    case 'bvp4c'
        if isfield(ocTrj.solverinfo,'yp')
            sol.solver=ocTrj.solver;
            sol.yp=ocTrj.solverinfo.yp;
        end
    case 'bvp5c'
        if isfield(ocTrj.solverinfo,'yp') &&  isfield(ocTrj.solverinfo,'ypmid')
            sol.x=ocTrj.x;
            sol.y=ocTrj.y;
            sol.solver=ocTrj.solver;
            sol.idata.yp=ocTrj.solverinfo.yp;
            sol.idata.ypmid=ocTrj.solverinfo.ypmid;
        end
    case 'bvp6c'
        if isfield(ocTrj.solverinfo,'yp') &&  isfield(ocTrj.solverinfo,'ypmid')
            sol.x=ocTrj.x;
            sol.y=ocTrj.y;
            sol.solver=ocTrj.solver;
            sol.yp=ocTrj.solverinfo.yp;
            sol.ypmid=ocTrj.solverinfo.ypmid;
        end
    case 'ode45'
        sol=ocTrj.solverinfo.ODESolveroutput;
    otherwise
        sol.solver='';
end
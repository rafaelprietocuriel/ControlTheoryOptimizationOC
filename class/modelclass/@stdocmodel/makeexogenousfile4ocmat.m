function modelfiles=makeexogenousfile4ocmat(ocObj,odes,moveflag)
%
% MAKEEXOGENOUSFILE4OCMAT generates files for exogenous dynamics:
% modelnameExogenousDynamics: containing the dynamics
% modelnameExogenousJacobian: containing the Jacobian of the dynamics with
% respect to the states, costates and exogenous state names
% modelnameExogenousParameterJacobian: containing the Jacobian of the
% dynamics with respect to the parameter names
%
% MAKEEXOGENOUSFILE4OCMAT(OCOBJ,ODES) standardmodel OCOBJ, symobolid ODES in the
% form sym('Dx=f(x)')
%
% MAKEEXOGENOUSFILE4OCMAT(OCOBJ,ODES,MOVEFLAG) MOVEFLAG=1 the files are moved to
% the model folder, MOVEFLAG=0 are stored in standard output folder. Default
% value is MOVEFLAG=1.
%
% Example usage:
% odes=sym('[DN=N*(b*l-delta),Dc=x^(sigma/beta)*h^(sigma*mu/beta+1)*(1-l-e)*N]')
% m=stdocmodel('skibajonesfinal');
% modelfiles=makeexogenousfile4ocmat(m,odes)

if nargin==2
    moveflag=1;
end
modelfiles=makeexogenousfile4ocmat(modelname(ocObj),odes,moveflag);
function J=equilibriumjacobian(odeObj,varargin)
%
% JACOBIAN returns the symbolic Jacobian or evaluated at a solution object.
% 
% J=JACOBIAN(ODEOBJ) ODEOBJ is a stdocmodel class. J is a string matrix of
% the formally written Jacobian for the canonical system.
% 
% J=JACOBIAN(ODEOBJ,SOLOBJ) SOLOBJ is a solution object
% (dynprimitive,octrajectory,ocasymptotic,occurve) a structure, including 
% at least the fieldnames 'y' and 'arcposition'. If SOLOBJ is empty the
% output is the same as for JACOBIAN(ODEOBJ). Otherwise the Jacobian
% evaluated at the solution object is returned. If SOLOBJ is an
% octrajectory consisting of multiple arcs J is a cell array of
% 3D-matrices, with the Jacobian values for each arc separately. 
%
% J=JACOBIAN(ODEOBJ,SOLOBJ,CONNECF) CONNECF=0 is the same as
% JACOBIAN(ODEOBJ,SOLOBJ). If CONNECF is 1 multiple arcs are returned
% "connected", i.e., the Jacobian values of all arcs are returned in one
% 3D-matrix.


J=[];
solstruct=[];
arcarg=[];
if isempty(odeObj)
    return
end
if nargin==1
    solstruct=[];
end
if nargin>=2
    solstruct=varargin{1};
end
if nargin>=3
    arcarg=varargin{2};
end
if isstruct(solstruct)
    try
        arcarg=solstruct.arcarg;
        indepvar=solstruct.x;
        solstruct=solstruct.y;
    catch
        ocmaterror('If the second input argument is a structure the fields ''y'', ''x'' and ''arcarg'' have to exist!')
    end
elseif isdynprimitive(solstruct) || isoctrajectory(solstruct) || isocasymptotic(solstruct) || isoccurve(solstruct)
    arcarg=arcargument(solstruct);
    %indepvar=time(odeObj,solstruct,1);
    indepvar=[];
    solstruct=dependentvar(solstruct);
end
if isempty(arcarg)
    arcarg=0;
end

if isempty(solstruct)
    % return symbolic expressions
    % if repflag=1 control values/Lagrange multipliers are replaced by the
    % corresponding symbolic expressions for the optimal values 
    J=feval(odeObj,'SymbolicEquilibriumEquationJacobian',arcarg);
else
    % return optimal control value evaluated at 'solstruct'
    %J=feval(odeObj,'EquilibriumEquationJacobian',0,solstruct,parametervalue(odeObj),arcarg);
    J=feval(odeObj,'EquilibriumEquationJacobian',solstruct,parametervalue(odeObj),arcarg);
end

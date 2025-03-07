function ocgEP = gradequilibrium(ocObj)
%
% EQUILIBRIUM returns the equilibria stored in OCOBJ
%
% OCEP=EQUILIBRIUM(OCOBJ) returns a cell arry of equilibria stored in the 
% 'Result' field 'Equilibrium' of OCOBJ.

% variable declaration
ocgEP=dynprimitive();

if isempty(ocObj)
    return
end

ocResult=result(ocObj);
if isfield(ocResult,'GradEquilibrium')
    ocgEP=ocResult.GradEquilibrium;
end
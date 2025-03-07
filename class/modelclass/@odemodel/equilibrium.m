function ocEP = equilibrium(odeObj)
%
% EQUILIBRIUM returns the equilibria stored in OCOBJ
%
% OCEP=EQUILIBRIUM(OCOBJ) returns a cell arry of equilibria stored in the 
% 'Result' field 'Equilibrium' of OCOBJ.

% variable declaration
ocEP=dynprimitive();

if isempty(odeObj)
    return
end

ocResult=result(odeObj);
if isfield(ocResult,'Equilibrium')
    ocEP=ocResult.Equilibrium;
end
function ocEP = equilibrium(dgObj)
%
% EQUILIBRIUM returns the equilibria stored in OCOBJ
%
% OCEP=EQUILIBRIUM(OCOBJ) returns a cell arry of equilibria stored in the 
% 'Result' field 'Equilibrium' of OCOBJ.

% variable declaration
ocEP=dynprimitive();

if isempty(dgObj)
    return
end

ocResult=result(dgObj);
if isfield(ocResult,'Equilibrium')
    ocEP=ocResult.Equilibrium;
end
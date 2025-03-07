function ocEP = fixpoint(ocObj)
%
% EQUILIBRIUM returns the equilibria stored in OCOBJ
%
% OCEP=EQUILIBRIUM(OCOBJ) returns a cell arry of equilibria stored in the 
% 'Result' field 'Equilibrium' of OCOBJ.

% variable declaration
ocEP=dynprimitive();

if isempty(ocObj)
    return
end

ocResult=result(ocObj);
if isfield(ocResult,'Fixpoint')
    ocEP=ocResult.Fixpoint;
end
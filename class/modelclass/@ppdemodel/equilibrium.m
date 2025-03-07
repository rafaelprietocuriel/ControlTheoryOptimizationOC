function ocEP = equilibrium(ppdeObj)
%
% EQUILIBRIUM returns the equilibria stored in ppdeObj
%
% OCEP=EQUILIBRIUM(ppdeObj) returns a cell arry of equilibria stored in the 
% 'Result' field 'Equilibrium' of ppdeObj.

% variable declaration
ocEP=ppdeprimitive();

if isempty(ppdeObj)
    return
end

ocResult=result(ppdeObj);
if isfield(ocResult,'Equilibrium')
    ocEP=ocResult.Equilibrium;
end
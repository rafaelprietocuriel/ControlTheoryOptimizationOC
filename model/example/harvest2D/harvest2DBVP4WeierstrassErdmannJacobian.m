function [Ja Jb Jpar]=harvest2DBVP4WeierstrassErdmannJacobian(dynVara,dynVarb,pararg,arcarg,switchid,switchcoord)
%
% the BC file for the Weierstrass Erdman conditions
                                                                                                 
% global variable
                                                                                                 
connecres=[];    % residuum for connecting hybrid paths
constrainres=[]; % residuum describing the boundary of constraints
multres=[];      % residuum describing the boundary of admissible Lagrange multipliers

switch switchid
    case ''
        Ja=[];
        Jb=[];
        Jpar=[];
        %ADDCASE
    otherwise
        error('switchid %s not defined',switchid)
end
out=[connecres;constrainres;multres];

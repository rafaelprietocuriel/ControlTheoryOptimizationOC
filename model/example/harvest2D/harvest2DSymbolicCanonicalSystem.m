function out = harvest2DSymbolicCanonicalSystem(arcarg)
%
% returns the symbolic canonical system of harvest2D model
                                                                                                           
u=harvest2DSymbolicOptimalControl(arcarg);
[lmmc lmsc]=harvest2DSymbolicLagrangeMultiplier(arcarg);
% Symbolic Canonical System
out=sym('[sigma*dynVar1*(1-dynVar1/m/dynVar2)-gamma/(C+tau)*dynVar1^theta1/(1+dynVar1^theta2)-eta*u1*dynVar1;(n-d*dynVar2-e*dynVar2*dynVar1)/epsilon;r*dynVar3-p*u1-dynVar3*(sigma*(1-dynVar1/m/dynVar2)-sigma*dynVar1/m/dynVar2-gamma/(C+tau)*dynVar1^theta1*theta1/dynVar1/(1+dynVar1^theta2)+gamma/(C+tau)*dynVar1^theta1/(1+dynVar1^theta2)^2*dynVar1^theta2*theta2/dynVar1-eta*u1)+dynVar4*e*dynVar2/epsilon;r*dynVar4-dynVar3*sigma*dynVar1^2/m/dynVar2^2-dynVar4*(-d-e*dynVar1)/epsilon]');
		out=subs(out,{'lmmc1'},{lmmc(1)});
%out=subs(out,{'u1'},{u(1)});

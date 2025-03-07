function J=cartel2dsymBCJacobian4Guard(depvara,depvarb,pararg,switchtime,arcid,edge,arc,arcoffset)
% returns the Jacobian of the guard boundary conditions for the saddle path continuation
%%
	% this file was automatically created: 24-Aug-2024 15:01:54
	% written by Dieter Grass, 2001 - 2024
	
if ~arcoffset
	J=-cartel2dsymBCJacobian4GuardArc(depvarb(:,arc),pararg,switchtime,edge(1,arc));
else
	J=cartel2dsymBCJacobian4GuardArc(depvara(:,arc+1),pararg,switchtime,edge(2,arc));
end
	
%--------------------------------------------------------------------------
function J=cartel2dsymBCJacobian4GuardArc(depvar,pararg,switchtime,arcid)
	
switch arcid
	case 0
		J=[0,0,0,0,0];
		J(1,1:4)=cartel2dsymDHamiltonianDx(switchtime,depvar,pararg,arcid);
	
	case 1
		J=[0,0,0,0,0];
		J(1,1:4)=cartel2dsymDHamiltonianDx(switchtime,depvar,pararg,arcid);
	
	case 2
		J=[0,0,0,0,0];
		J(1,1:4)=cartel2dsymDHamiltonianDx(switchtime,depvar,pararg,arcid);
	
	case 3
		J=[0,0,0,0];
		J(1,1:4)=cartel2dsymDHamiltonianDx(switchtime,depvar,pararg,arcid);
	
	case 4
		J=[0,0,0,0,0,0];
		J(1,1:4)=cartel2dsymDHamiltonianDx(switchtime,depvar,pararg,arcid);
	
	case 5
		J=[0,0,0,0,0];
		J(1,1:4)=cartel2dsymDHamiltonianDx(switchtime,depvar,pararg,arcid);
	
	case 6
		J=[0,0,0,0];
		J(1,1:4)=cartel2dsymDHamiltonianDx(switchtime,depvar,pararg,arcid);
	
	case 7
		J=[0,0,0,0,0];
		J(1,1:4)=cartel2dsymDHamiltonianDx(switchtime,depvar,pararg,arcid);
	
	case 8
		J=[0,0,0,0];
		J(1,1:4)=cartel2dsymDHamiltonianDx(switchtime,depvar,pararg,arcid);
	
	case 9
		J=[0,0,0,0,0];
		J(1,1:4)=cartel2dsymDHamiltonianDx(switchtime,depvar,pararg,arcid);
	
	case 10
		J=[0,0,0,0,0];
		J(1,1:4)=cartel2dsymDHamiltonianDx(switchtime,depvar,pararg,arcid);
	
	case 11
		J=[0,0,0,0];
		J(1,1:4)=cartel2dsymDHamiltonianDx(switchtime,depvar,pararg,arcid);
	
	case 12
		J=[0,0,0,0];
		J(1,1:4)=cartel2dsymDHamiltonianDx(switchtime,depvar,pararg,arcid);
	
	case 13
		J=[0,0,0,0];
		J(1,1:4)=cartel2dsymDHamiltonianDx(switchtime,depvar,pararg,arcid);
	
	case 14
		J=[0,0,0,0];
		J(1,1:4)=cartel2dsymDHamiltonianDx(switchtime,depvar,pararg,arcid);
	
end

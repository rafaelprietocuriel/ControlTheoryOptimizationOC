function J=cartel2dsymBCJacobian4Reset(depvara,depvarb,pararg,freepar,arcid,edge,arc,arcoffset)
% returns the Jacobian of the guard boundary conditions for the saddle path continuation
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
	
if ~arcoffset
	J=-cartel2dsymBCJacobian4ResetArc(depvarb(:,arc),pararg,freepar,edge(1));
else
	J=cartel2dsymBCJacobian4ResetArc(depvara(:,arc),pararg,freepar,edge(2));
end
	
%--------------------------------------------------------------------------
function J=cartel2dsymBCJacobian4ResetArc(depvar,pararg,freepar,arcid)
	
switch arcid
	case 0
		J=zeros(4,5);
		J(1:4,1:4)=eye(4);
	
	case 1
		J=zeros(4,5);
		J(1:4,1:4)=eye(4);
	
	case 2
		J=zeros(4,5);
		J(1:4,1:4)=eye(4);
	
	case 3
		J=zeros(4,4);
		J(1:4,1:4)=eye(4);
	
	case 4
		J=zeros(4,6);
		J(1:4,1:4)=eye(4);
	
	case 5
		J=zeros(4,5);
		J(1:4,1:4)=eye(4);
	
	case 6
		J=zeros(4,4);
		J(1:4,1:4)=eye(4);
	
	case 7
		J=zeros(4,5);
		J(1:4,1:4)=eye(4);
	
	case 8
		J=zeros(4,4);
		J(1:4,1:4)=eye(4);
	
	case 9
		J=zeros(4,5);
		J(1:4,1:4)=eye(4);
	
	case 10
		J=zeros(4,5);
		J(1:4,1:4)=eye(4);
	
	case 11
		J=zeros(4,4);
		J(1:4,1:4)=eye(4);
	
	case 12
		J=zeros(4,4);
		J(1:4,1:4)=eye(4);
	
	case 13
		J=zeros(4,4);
		J(1:4,1:4)=eye(4);
	
	case 14
		J=zeros(4,4);
		J(1:4,1:4)=eye(4);
	
end

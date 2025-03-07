function out=cartel2dsymObjectiveFunctionDerivativeTime(t,depvar,par,arcid)
% returns the time derivative of the discounted objective function
%
% this file was automatically created: 24-Aug-2024 15:01:54
% written by Dieter Grass, 2001 - 2024
	
r=par(1);
beta=par(2);
alpha=par(3);
rho=par(4);
theta=par(5);
omega=par(6);
delta=par(7);
h=par(8);
gamma1=par(9);
gamma2=par(10);
kappa=par(11);
eta=par(12);
ulow=par(13);
vlow=par(14);
w=par(15);
s=par(16);
b=par(17);
xi=par(18);
epsilon=par(19);
tau=par(20);
B=par(21);
	
switch arcid
	case 0
		out=-r*exp(-r*t)*(-w*s*depvar(1,:)*depvar(2,:)-(1-w)*h*(depvar(1,:)+depvar(2,:))-gamma1*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))-gamma2*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))-kappa*depvar(5,:));
	
	case 1
		out=-r*exp(-r*t)*(-w*s*depvar(1,:)*depvar(2,:)-(1-w)*h*(depvar(1,:)+depvar(2,:))-gamma1*ulow-gamma2*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))-kappa*depvar(5,:));
	
	case 2
		out=-r*exp(-r*t)*(-w*s*depvar(1,:)*depvar(2,:)-(1-w)*h*(depvar(1,:)+depvar(2,:))-gamma1*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))-gamma2*ulow-kappa*depvar(5,:));
	
	case 3
		out=-r*exp(-r*t)*(-w*s*depvar(1,:)*depvar(2,:)-(1-w)*h*(depvar(1,:)+depvar(2,:))-gamma1*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))-gamma2*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))-kappa*vlow);
	
	case 4
		out=-r*exp(-r*t)*(-w*s*depvar(1,:)*depvar(2,:)-(1-w)*h*(depvar(1,:)+depvar(2,:))-gamma1*depvar(5,:)-gamma2*depvar(6,:)-kappa*(B-depvar(5,:)-depvar(6,:)));
	
	case 5
		out=-r*exp(-r*t)*(-w*s*depvar(1,:)*depvar(2,:)-(1-w)*h*(depvar(1,:)+depvar(2,:))-gamma1*ulow-gamma2*ulow-kappa*depvar(5,:));
	
	case 6
		out=-r*exp(-r*t)*(-w*s*depvar(1,:)*depvar(2,:)-(1-w)*h*(depvar(1,:)+depvar(2,:))-gamma1*ulow-gamma2*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))-kappa*vlow);
	
	case 7
		out=-r*exp(-r*t)*(-w*s*depvar(1,:)*depvar(2,:)-(1-w)*h*(depvar(1,:)+depvar(2,:))-gamma1*ulow-gamma2*depvar(5,:)-kappa*(B-ulow-depvar(5,:)));
	
	case 8
		out=-r*exp(-r*t)*(-w*s*depvar(1,:)*depvar(2,:)-(1-w)*h*(depvar(1,:)+depvar(2,:))-gamma1*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))-gamma2*ulow-kappa*vlow);
	
	case 9
		out=-r*exp(-r*t)*(-w*s*depvar(1,:)*depvar(2,:)-(1-w)*h*(depvar(1,:)+depvar(2,:))-gamma1*depvar(5,:)-gamma2*ulow-kappa*(B-depvar(5,:)-ulow));
	
	case 10
		out=-r*exp(-r*t)*(-w*s*depvar(1,:)*depvar(2,:)-(1-w)*h*(depvar(1,:)+depvar(2,:))-gamma1*(B-depvar(5,:)-vlow)-gamma2*depvar(5,:)-kappa*vlow);
	
	case 11
		out=-r*exp(-r*t)*(-w*s*depvar(1,:)*depvar(2,:)-(1-w)*h*(depvar(1,:)+depvar(2,:))-gamma1*ulow-gamma2*ulow-kappa*vlow);
	
	case 12
		out=-r*exp(-r*t)*(-w*s*depvar(1,:)*depvar(2,:)-(1-w)*h*(depvar(1,:)+depvar(2,:))-gamma1*ulow-gamma2*ulow-kappa*(B-2*ulow));
	
	case 13
		out=-r*exp(-r*t)*(-w*s*depvar(1,:)*depvar(2,:)-(1-w)*h*(depvar(1,:)+depvar(2,:))-gamma1*ulow-gamma2*(B-ulow-vlow)-kappa*vlow);
	
	case 14
		out=-r*exp(-r*t)*(-w*s*depvar(1,:)*depvar(2,:)-(1-w)*h*(depvar(1,:)+depvar(2,:))-gamma1*(B-ulow-vlow)-gamma2*ulow-kappa*vlow);
	
end

function out=cartel2dsymDHamiltonianDx(t,depvar,par,arcid)
% returns the derivative of the Hamiltonian with respect to the (co)state variables cartel2dsym
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
	
ctrl=cartel2dsymOptimalControl(t,depvar,par,arcid);
[lagmcc,lagmsc]=cartel2dsymLagrangeMultiplier(t,depvar,par,arcid);
	
switch arcid
	case 0
		out=[-w*s*depvar(2,:)-(1-w)*h+gamma1*beta/(depvar(1,:)+epsilon)/(-1+alpha)*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))+depvar(3,:)*((b+(delta-b)*exp(-xi*depvar(5,:)))*rho+eta*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))^alpha*alpha*beta/(depvar(1,:)+epsilon)/(-1+alpha)*(depvar(1,:)+epsilon)^beta-eta*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))^alpha*(depvar(1,:)+epsilon)^beta*beta/(depvar(1,:)+epsilon)-theta*depvar(2,:)-2*omega*depvar(1,:))-depvar(4,:)*theta*depvar(2,:), ...
			 -w*s*depvar(1,:)-(1-w)*h+gamma2*beta/(depvar(2,:)+epsilon)/(-1+alpha)*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))-depvar(3,:)*theta*depvar(1,:)+depvar(4,:)*((b+(delta-b)*exp(-xi*depvar(5,:)))*rho+eta*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))^alpha*alpha*beta/(depvar(2,:)+epsilon)/(-1+alpha)*(depvar(2,:)+epsilon)^beta-eta*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))^alpha*(depvar(2,:)+epsilon)^beta*beta/(depvar(2,:)+epsilon)-theta*depvar(1,:)-2*omega*depvar(2,:)), ...
			 gamma1/depvar(3,:)/(-1+alpha)*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))+tau+(b+(delta-b)*exp(-xi*depvar(5,:)))*rho*depvar(1,:)-eta*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))^alpha*(depvar(1,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(1,:)^2+eta*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))^alpha*alpha/(-1+alpha)*(depvar(1,:)+epsilon)^beta, ...
			 gamma2/depvar(4,:)/(-1+alpha)*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))+tau+(b+(delta-b)*exp(-xi*depvar(5,:)))*rho*depvar(2,:)-eta*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))^alpha*(depvar(2,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(2,:)^2+eta*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))^alpha*alpha/(-1+alpha)*(depvar(2,:)+epsilon)^beta];
	
	case 1
		out=[-w*s*depvar(2,:)-(1-w)*h+depvar(3,:)*((b+(delta-b)*exp(-xi*depvar(5,:)))*rho-eta*ulow^alpha*(depvar(1,:)+epsilon)^beta*beta/(depvar(1,:)+epsilon)-theta*depvar(2,:)-2*omega*depvar(1,:))-depvar(4,:)*theta*depvar(2,:), ...
			 -w*s*depvar(1,:)-(1-w)*h+gamma2*beta/(depvar(2,:)+epsilon)/(-1+alpha)*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))-depvar(3,:)*theta*depvar(1,:)+depvar(4,:)*((b+(delta-b)*exp(-xi*depvar(5,:)))*rho+eta*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))^alpha*alpha*beta/(depvar(2,:)+epsilon)/(-1+alpha)*(depvar(2,:)+epsilon)^beta-eta*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))^alpha*(depvar(2,:)+epsilon)^beta*beta/(depvar(2,:)+epsilon)-theta*depvar(1,:)-2*omega*depvar(2,:)), ...
			 tau+(b+(delta-b)*exp(-xi*depvar(5,:)))*rho*depvar(1,:)-eta*ulow^alpha*(depvar(1,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(1,:)^2, ...
			 gamma2/depvar(4,:)/(-1+alpha)*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))+tau+(b+(delta-b)*exp(-xi*depvar(5,:)))*rho*depvar(2,:)-eta*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))^alpha*(depvar(2,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(2,:)^2+eta*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))^alpha*alpha/(-1+alpha)*(depvar(2,:)+epsilon)^beta];
	
	case 2
		out=[-w*s*depvar(2,:)-(1-w)*h+gamma1*beta/(depvar(1,:)+epsilon)/(-1+alpha)*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))+depvar(3,:)*((b+(delta-b)*exp(-xi*depvar(5,:)))*rho+eta*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))^alpha*alpha*beta/(depvar(1,:)+epsilon)/(-1+alpha)*(depvar(1,:)+epsilon)^beta-eta*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))^alpha*(depvar(1,:)+epsilon)^beta*beta/(depvar(1,:)+epsilon)-theta*depvar(2,:)-2*omega*depvar(1,:))-depvar(4,:)*theta*depvar(2,:), ...
			 -w*s*depvar(1,:)-(1-w)*h-depvar(3,:)*theta*depvar(1,:)+depvar(4,:)*((b+(delta-b)*exp(-xi*depvar(5,:)))*rho-eta*ulow^alpha*(depvar(2,:)+epsilon)^beta*beta/(depvar(2,:)+epsilon)-theta*depvar(1,:)-2*omega*depvar(2,:)), ...
			 gamma1/depvar(3,:)/(-1+alpha)*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))+tau+(b+(delta-b)*exp(-xi*depvar(5,:)))*rho*depvar(1,:)-eta*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))^alpha*(depvar(1,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(1,:)^2+eta*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))^alpha*alpha/(-1+alpha)*(depvar(1,:)+epsilon)^beta, ...
			 tau+(b+(delta-b)*exp(-xi*depvar(5,:)))*rho*depvar(2,:)-eta*ulow^alpha*(depvar(2,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(2,:)^2];
	
	case 3
		out=[-w*s*depvar(2,:)-(1-w)*h+gamma1*beta/(depvar(1,:)+epsilon)/(-1+alpha)*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))+depvar(3,:)*((b+(delta-b)*exp(-xi*vlow))*rho+eta*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))^alpha*alpha*beta/(depvar(1,:)+epsilon)/(-1+alpha)*(depvar(1,:)+epsilon)^beta-eta*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))^alpha*(depvar(1,:)+epsilon)^beta*beta/(depvar(1,:)+epsilon)-theta*depvar(2,:)-2*omega*depvar(1,:))-depvar(4,:)*theta*depvar(2,:), ...
			 -w*s*depvar(1,:)-(1-w)*h+gamma2*beta/(depvar(2,:)+epsilon)/(-1+alpha)*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))-depvar(3,:)*theta*depvar(1,:)+depvar(4,:)*((b+(delta-b)*exp(-xi*vlow))*rho+eta*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))^alpha*alpha*beta/(depvar(2,:)+epsilon)/(-1+alpha)*(depvar(2,:)+epsilon)^beta-eta*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))^alpha*(depvar(2,:)+epsilon)^beta*beta/(depvar(2,:)+epsilon)-theta*depvar(1,:)-2*omega*depvar(2,:)), ...
			 gamma1/depvar(3,:)/(-1+alpha)*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))+tau+(b+(delta-b)*exp(-xi*vlow))*rho*depvar(1,:)-eta*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))^alpha*(depvar(1,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(1,:)^2+eta*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))^alpha*alpha/(-1+alpha)*(depvar(1,:)+epsilon)^beta, ...
			 gamma2/depvar(4,:)/(-1+alpha)*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))+tau+(b+(delta-b)*exp(-xi*vlow))*rho*depvar(2,:)-eta*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))^alpha*(depvar(2,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(2,:)^2+eta*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))^alpha*alpha/(-1+alpha)*(depvar(2,:)+epsilon)^beta];
	
	case 4
		out=[-w*s*depvar(2,:)-(1-w)*h+depvar(3,:)*((b+(delta-b)*exp(-xi*(B-depvar(5,:)-depvar(6,:))))*rho-eta*depvar(5,:)^alpha*(depvar(1,:)+epsilon)^beta*beta/(depvar(1,:)+epsilon)-theta*depvar(2,:)-2*omega*depvar(1,:))-depvar(4,:)*theta*depvar(2,:), ...
			 -w*s*depvar(1,:)-(1-w)*h-depvar(3,:)*theta*depvar(1,:)+depvar(4,:)*((b+(delta-b)*exp(-xi*(B-depvar(5,:)-depvar(6,:))))*rho-eta*depvar(6,:)^alpha*(depvar(2,:)+epsilon)^beta*beta/(depvar(2,:)+epsilon)-theta*depvar(1,:)-2*omega*depvar(2,:)), ...
			 tau+(b+(delta-b)*exp(-xi*(B-depvar(5,:)-depvar(6,:))))*rho*depvar(1,:)-eta*depvar(5,:)^alpha*(depvar(1,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(1,:)^2, ...
			 tau+(b+(delta-b)*exp(-xi*(B-depvar(5,:)-depvar(6,:))))*rho*depvar(2,:)-eta*depvar(6,:)^alpha*(depvar(2,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(2,:)^2];
	
	case 5
		out=[-w*s*depvar(2,:)-(1-w)*h+depvar(3,:)*((b+(delta-b)*exp(-xi*depvar(5,:)))*rho-eta*ulow^alpha*(depvar(1,:)+epsilon)^beta*beta/(depvar(1,:)+epsilon)-theta*depvar(2,:)-2*omega*depvar(1,:))-depvar(4,:)*theta*depvar(2,:), ...
			 -w*s*depvar(1,:)-(1-w)*h-depvar(3,:)*theta*depvar(1,:)+depvar(4,:)*((b+(delta-b)*exp(-xi*depvar(5,:)))*rho-eta*ulow^alpha*(depvar(2,:)+epsilon)^beta*beta/(depvar(2,:)+epsilon)-theta*depvar(1,:)-2*omega*depvar(2,:)), ...
			 tau+(b+(delta-b)*exp(-xi*depvar(5,:)))*rho*depvar(1,:)-eta*ulow^alpha*(depvar(1,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(1,:)^2, ...
			 tau+(b+(delta-b)*exp(-xi*depvar(5,:)))*rho*depvar(2,:)-eta*ulow^alpha*(depvar(2,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(2,:)^2];
	
	case 6
		out=[-w*s*depvar(2,:)-(1-w)*h+depvar(3,:)*((b+(delta-b)*exp(-xi*vlow))*rho-eta*ulow^alpha*(depvar(1,:)+epsilon)^beta*beta/(depvar(1,:)+epsilon)-theta*depvar(2,:)-2*omega*depvar(1,:))-depvar(4,:)*theta*depvar(2,:), ...
			 -w*s*depvar(1,:)-(1-w)*h+gamma2*beta/(depvar(2,:)+epsilon)/(-1+alpha)*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))-depvar(3,:)*theta*depvar(1,:)+depvar(4,:)*((b+(delta-b)*exp(-xi*vlow))*rho+eta*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))^alpha*alpha*beta/(depvar(2,:)+epsilon)/(-1+alpha)*(depvar(2,:)+epsilon)^beta-eta*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))^alpha*(depvar(2,:)+epsilon)^beta*beta/(depvar(2,:)+epsilon)-theta*depvar(1,:)-2*omega*depvar(2,:)), ...
			 tau+(b+(delta-b)*exp(-xi*vlow))*rho*depvar(1,:)-eta*ulow^alpha*(depvar(1,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(1,:)^2, ...
			 gamma2/depvar(4,:)/(-1+alpha)*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))+tau+(b+(delta-b)*exp(-xi*vlow))*rho*depvar(2,:)-eta*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))^alpha*(depvar(2,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(2,:)^2+eta*exp((log(-gamma2/depvar(4,:)/eta/alpha)-beta*log(depvar(2,:)+epsilon))/(-1+alpha))^alpha*alpha/(-1+alpha)*(depvar(2,:)+epsilon)^beta];
	
	case 7
		out=[-w*s*depvar(2,:)-(1-w)*h+depvar(3,:)*((b+(delta-b)*exp(-xi*(B-ulow-depvar(5,:))))*rho-eta*ulow^alpha*(depvar(1,:)+epsilon)^beta*beta/(depvar(1,:)+epsilon)-theta*depvar(2,:)-2*omega*depvar(1,:))-depvar(4,:)*theta*depvar(2,:), ...
			 -w*s*depvar(1,:)-(1-w)*h-depvar(3,:)*theta*depvar(1,:)+depvar(4,:)*((b+(delta-b)*exp(-xi*(B-ulow-depvar(5,:))))*rho-eta*depvar(5,:)^alpha*(depvar(2,:)+epsilon)^beta*beta/(depvar(2,:)+epsilon)-theta*depvar(1,:)-2*omega*depvar(2,:)), ...
			 tau+(b+(delta-b)*exp(-xi*(B-ulow-depvar(5,:))))*rho*depvar(1,:)-eta*ulow^alpha*(depvar(1,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(1,:)^2, ...
			 tau+(b+(delta-b)*exp(-xi*(B-ulow-depvar(5,:))))*rho*depvar(2,:)-eta*depvar(5,:)^alpha*(depvar(2,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(2,:)^2];
	
	case 8
		out=[-w*s*depvar(2,:)-(1-w)*h+gamma1*beta/(depvar(1,:)+epsilon)/(-1+alpha)*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))+depvar(3,:)*((b+(delta-b)*exp(-xi*vlow))*rho+eta*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))^alpha*alpha*beta/(depvar(1,:)+epsilon)/(-1+alpha)*(depvar(1,:)+epsilon)^beta-eta*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))^alpha*(depvar(1,:)+epsilon)^beta*beta/(depvar(1,:)+epsilon)-theta*depvar(2,:)-2*omega*depvar(1,:))-depvar(4,:)*theta*depvar(2,:), ...
			 -w*s*depvar(1,:)-(1-w)*h-depvar(3,:)*theta*depvar(1,:)+depvar(4,:)*((b+(delta-b)*exp(-xi*vlow))*rho-eta*ulow^alpha*(depvar(2,:)+epsilon)^beta*beta/(depvar(2,:)+epsilon)-theta*depvar(1,:)-2*omega*depvar(2,:)), ...
			 gamma1/depvar(3,:)/(-1+alpha)*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))+tau+(b+(delta-b)*exp(-xi*vlow))*rho*depvar(1,:)-eta*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))^alpha*(depvar(1,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(1,:)^2+eta*exp(-(-log(-gamma1/depvar(3,:)/eta/alpha)+beta*log(depvar(1,:)+epsilon))/(-1+alpha))^alpha*alpha/(-1+alpha)*(depvar(1,:)+epsilon)^beta, ...
			 tau+(b+(delta-b)*exp(-xi*vlow))*rho*depvar(2,:)-eta*ulow^alpha*(depvar(2,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(2,:)^2];
	
	case 9
		out=[-w*s*depvar(2,:)-(1-w)*h+depvar(3,:)*((b+(delta-b)*exp(-xi*(B-depvar(5,:)-ulow)))*rho-eta*depvar(5,:)^alpha*(depvar(1,:)+epsilon)^beta*beta/(depvar(1,:)+epsilon)-theta*depvar(2,:)-2*omega*depvar(1,:))-depvar(4,:)*theta*depvar(2,:), ...
			 -w*s*depvar(1,:)-(1-w)*h-depvar(3,:)*theta*depvar(1,:)+depvar(4,:)*((b+(delta-b)*exp(-xi*(B-depvar(5,:)-ulow)))*rho-eta*ulow^alpha*(depvar(2,:)+epsilon)^beta*beta/(depvar(2,:)+epsilon)-theta*depvar(1,:)-2*omega*depvar(2,:)), ...
			 tau+(b+(delta-b)*exp(-xi*(B-depvar(5,:)-ulow)))*rho*depvar(1,:)-eta*depvar(5,:)^alpha*(depvar(1,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(1,:)^2, ...
			 tau+(b+(delta-b)*exp(-xi*(B-depvar(5,:)-ulow)))*rho*depvar(2,:)-eta*ulow^alpha*(depvar(2,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(2,:)^2];
	
	case 10
		out=[-w*s*depvar(2,:)-(1-w)*h+depvar(3,:)*((b+(delta-b)*exp(-xi*vlow))*rho-eta*(B-depvar(5,:)-vlow)^alpha*(depvar(1,:)+epsilon)^beta*beta/(depvar(1,:)+epsilon)-theta*depvar(2,:)-2*omega*depvar(1,:))-depvar(4,:)*theta*depvar(2,:), ...
			 -w*s*depvar(1,:)-(1-w)*h-depvar(3,:)*theta*depvar(1,:)+depvar(4,:)*((b+(delta-b)*exp(-xi*vlow))*rho-eta*depvar(5,:)^alpha*(depvar(2,:)+epsilon)^beta*beta/(depvar(2,:)+epsilon)-theta*depvar(1,:)-2*omega*depvar(2,:)), ...
			 tau+(b+(delta-b)*exp(-xi*vlow))*rho*depvar(1,:)-eta*(B-depvar(5,:)-vlow)^alpha*(depvar(1,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(1,:)^2, ...
			 tau+(b+(delta-b)*exp(-xi*vlow))*rho*depvar(2,:)-eta*depvar(5,:)^alpha*(depvar(2,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(2,:)^2];
	
	case 11
		out=[-w*s*depvar(2,:)-(1-w)*h+depvar(3,:)*((b+(delta-b)*exp(-xi*vlow))*rho-eta*ulow^alpha*(depvar(1,:)+epsilon)^beta*beta/(depvar(1,:)+epsilon)-theta*depvar(2,:)-2*omega*depvar(1,:))-depvar(4,:)*theta*depvar(2,:), ...
			 -w*s*depvar(1,:)-(1-w)*h-depvar(3,:)*theta*depvar(1,:)+depvar(4,:)*((b+(delta-b)*exp(-xi*vlow))*rho-eta*ulow^alpha*(depvar(2,:)+epsilon)^beta*beta/(depvar(2,:)+epsilon)-theta*depvar(1,:)-2*omega*depvar(2,:)), ...
			 tau+(b+(delta-b)*exp(-xi*vlow))*rho*depvar(1,:)-eta*ulow^alpha*(depvar(1,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(1,:)^2, ...
			 tau+(b+(delta-b)*exp(-xi*vlow))*rho*depvar(2,:)-eta*ulow^alpha*(depvar(2,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(2,:)^2];
	
	case 12
		out=[-w*s*depvar(2,:)-(1-w)*h+depvar(3,:)*((b+(delta-b)*exp(-xi*(B-2*ulow)))*rho-eta*ulow^alpha*(depvar(1,:)+epsilon)^beta*beta/(depvar(1,:)+epsilon)-theta*depvar(2,:)-2*omega*depvar(1,:))-depvar(4,:)*theta*depvar(2,:), ...
			 -w*s*depvar(1,:)-(1-w)*h-depvar(3,:)*theta*depvar(1,:)+depvar(4,:)*((b+(delta-b)*exp(-xi*(B-2*ulow)))*rho-eta*ulow^alpha*(depvar(2,:)+epsilon)^beta*beta/(depvar(2,:)+epsilon)-theta*depvar(1,:)-2*omega*depvar(2,:)), ...
			 tau+(b+(delta-b)*exp(-xi*(B-2*ulow)))*rho*depvar(1,:)-eta*ulow^alpha*(depvar(1,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(1,:)^2, ...
			 tau+(b+(delta-b)*exp(-xi*(B-2*ulow)))*rho*depvar(2,:)-eta*ulow^alpha*(depvar(2,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(2,:)^2];
	
	case 13
		out=[-w*s*depvar(2,:)-(1-w)*h+depvar(3,:)*((b+(delta-b)*exp(-xi*vlow))*rho-eta*ulow^alpha*(depvar(1,:)+epsilon)^beta*beta/(depvar(1,:)+epsilon)-theta*depvar(2,:)-2*omega*depvar(1,:))-depvar(4,:)*theta*depvar(2,:), ...
			 -w*s*depvar(1,:)-(1-w)*h-depvar(3,:)*theta*depvar(1,:)+depvar(4,:)*((b+(delta-b)*exp(-xi*vlow))*rho-eta*(B-ulow-vlow)^alpha*(depvar(2,:)+epsilon)^beta*beta/(depvar(2,:)+epsilon)-theta*depvar(1,:)-2*omega*depvar(2,:)), ...
			 tau+(b+(delta-b)*exp(-xi*vlow))*rho*depvar(1,:)-eta*ulow^alpha*(depvar(1,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(1,:)^2, ...
			 tau+(b+(delta-b)*exp(-xi*vlow))*rho*depvar(2,:)-eta*(B-ulow-vlow)^alpha*(depvar(2,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(2,:)^2];
	
	case 14
		out=[-w*s*depvar(2,:)-(1-w)*h+depvar(3,:)*((b+(delta-b)*exp(-xi*vlow))*rho-eta*(B-ulow-vlow)^alpha*(depvar(1,:)+epsilon)^beta*beta/(depvar(1,:)+epsilon)-theta*depvar(2,:)-2*omega*depvar(1,:))-depvar(4,:)*theta*depvar(2,:), ...
			 -w*s*depvar(1,:)-(1-w)*h-depvar(3,:)*theta*depvar(1,:)+depvar(4,:)*((b+(delta-b)*exp(-xi*vlow))*rho-eta*ulow^alpha*(depvar(2,:)+epsilon)^beta*beta/(depvar(2,:)+epsilon)-theta*depvar(1,:)-2*omega*depvar(2,:)), ...
			 tau+(b+(delta-b)*exp(-xi*vlow))*rho*depvar(1,:)-eta*(B-ulow-vlow)^alpha*(depvar(1,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(1,:)^2, ...
			 tau+(b+(delta-b)*exp(-xi*vlow))*rho*depvar(2,:)-eta*ulow^alpha*(depvar(2,:)+epsilon)^beta-theta*depvar(1,:)*depvar(2,:)-omega*depvar(2,:)^2];
	
end

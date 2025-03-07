function out = harvest2DHamiltonian(t,dynVar,pararg,arcarg)
%
% returns the Hamiltonian of harvest2D model
                                                                                                                 
% Parameter Value
r=pararg(1);      
C=pararg(2);      
gamma=pararg(3);  
d=pararg(4);      
e=pararg(5);      
epsilon=pararg(6);
eta=pararg(7);    
m=pararg(8);      
n=pararg(9);      
p=pararg(10);     
phi=pararg(11);   
sigma=pararg(12); 
tau=pararg(13);   
theta1=pararg(14);
theta2=pararg(15);
ulow=pararg(16);  
T=pararg(17);     
% Functions

% Hamiltonian
switch arcarg
    case 1
        u=dynVar(5,:);
    case 2
        u=harvest2DOptimalControl(t,dynVar,pararg,arcarg);
end
out=[p.*u(1,:).*dynVar(1,:)-phi.*u(1,:).^2+dynVar(3,:).*(sigma.*dynVar(1,:).*(1-dynVar(1,:)./m./dynVar(2,:))-gamma./(C+tau).*dynVar(1,:).^theta1./(1+dynVar(1,:).^theta2)-eta.*u(1,:).*dynVar(1,:))+dynVar(4,:).*(n-d.*dynVar(2,:)-e.*dynVar(2,:).*dynVar(1,:))./epsilon];

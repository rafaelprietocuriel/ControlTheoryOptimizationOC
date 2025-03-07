function out = harvest2DDHamiltonianDu(t,dynVar,pararg,arcarg)
%
% returns the partial derivative of the Hamiltonian with repsect to
% the control variables
                                                                                            
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
                                                                                            
u=dynVar(5,:);
[lmmc lmsc]=harvest2DLagrangeMultiplier(t,dynVar,pararg,arcarg);
out=[p.*dynVar(1,:)-2.*phi.*u(1,:)-dynVar(3,:).*eta.*dynVar(1,:)+lmmc(1,:)];

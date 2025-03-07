function J=harvest2DCanonicalSystemICJacobian(t,dynVar,pararg,arcarg)
%
% the dynamics file for model harvest2D for the BVP
                                                                                                                                        
% Jacobian system

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
                                                                                            
% Canonical System
switch arcarg
    case 1
        u=dynVar(5);
        J=[sigma*(1-dynVar(1)/m/dynVar(2))-sigma*dynVar(1)/m/dynVar(2)-gamma/(C+tau)*dynVar(1)^theta1*theta1/dynVar(1)/(1+dynVar(1)^theta2)+gamma/(C+tau)*dynVar(1)^theta1/(1+dynVar(1)^theta2)^2*dynVar(1)^theta2*theta2/dynVar(1)-eta*u(1), sigma*dynVar(1)^2/m/dynVar(2)^2, 0,0,-eta*dynVar(1); ...
            -e*dynVar(2)/epsilon, (-d-e*dynVar(1))/epsilon,0 ,0 ,0; ...
            -dynVar(3)*(-2*sigma/m/dynVar(2)-gamma/(C+tau)*dynVar(1)^theta1*theta1^2/dynVar(1)^2/(1+dynVar(1)^theta2)+gamma/(C+tau)*dynVar(1)^theta1*theta1/dynVar(1)^2/(1+dynVar(1)^theta2)+2*gamma/(C+tau)*dynVar(1)^theta1*theta1/dynVar(1)^2/(1+dynVar(1)^theta2)^2*dynVar(1)^theta2*theta2-2*gamma/(C+tau)*dynVar(1)^theta1/(1+dynVar(1)^theta2)^3*(dynVar(1)^theta2)^2*theta2^2/dynVar(1)^2+gamma/(C+tau)*dynVar(1)^theta1/(1+dynVar(1)^theta2)^2*dynVar(1)^theta2*theta2^2/dynVar(1)^2-gamma/(C+tau)*dynVar(1)^theta1/(1+dynVar(1)^theta2)^2*dynVar(1)^theta2*theta2/dynVar(1)^2), -2*dynVar(3)*sigma*dynVar(1)/m/dynVar(2)^2+dynVar(4)*e/epsilon, r-sigma*(1-dynVar(1)/m/dynVar(2))+sigma*dynVar(1)/m/dynVar(2)+gamma/(C+tau)*dynVar(1)^theta1*theta1/dynVar(1)/(1+dynVar(1)^theta2)-gamma/(C+tau)*dynVar(1)^theta1/(1+dynVar(1)^theta2)^2*dynVar(1)^theta2*theta2/dynVar(1)+eta*u(1), e*dynVar(2)/epsilon, -p+dynVar(3)*eta; ...
            -2*dynVar(3)*sigma*dynVar(1)/m/dynVar(2)^2+dynVar(4)*e/epsilon, 2*dynVar(3)*sigma*dynVar(1)^2/m/dynVar(2)^3, -sigma*dynVar(1)^2/m/dynVar(2)^2, r-(-d-e*dynVar(1))/epsilon, 0; ...
            p-dynVar(3)*eta, 0, -eta*dynVar(1), 0, -2*phi];
    case 2
        u=harvest2DOptimalControl(t,dynVar,pararg,arcarg);
        J=[sigma*(1-dynVar(1)/m/dynVar(2))-sigma*dynVar(1)/m/dynVar(2)-gamma/(C+tau)*dynVar(1)^theta1*theta1/dynVar(1)/(1+dynVar(1)^theta2)+gamma/(C+tau)*dynVar(1)^theta1/(1+dynVar(1)^theta2)^2*dynVar(1)^theta2*theta2/dynVar(1)-eta*u(1), sigma*dynVar(1)^2/m/dynVar(2)^2, 0,0; ...
            -e*dynVar(2)/epsilon, (-d-e*dynVar(1))/epsilon,0 ,0; ...
            -dynVar(3)*(-2*sigma/m/dynVar(2)-gamma/(C+tau)*dynVar(1)^theta1*theta1^2/dynVar(1)^2/(1+dynVar(1)^theta2)+gamma/(C+tau)*dynVar(1)^theta1*theta1/dynVar(1)^2/(1+dynVar(1)^theta2)+2*gamma/(C+tau)*dynVar(1)^theta1*theta1/dynVar(1)^2/(1+dynVar(1)^theta2)^2*dynVar(1)^theta2*theta2-2*gamma/(C+tau)*dynVar(1)^theta1/(1+dynVar(1)^theta2)^3*(dynVar(1)^theta2)^2*theta2^2/dynVar(1)^2+gamma/(C+tau)*dynVar(1)^theta1/(1+dynVar(1)^theta2)^2*dynVar(1)^theta2*theta2^2/dynVar(1)^2-gamma/(C+tau)*dynVar(1)^theta1/(1+dynVar(1)^theta2)^2*dynVar(1)^theta2*theta2/dynVar(1)^2), -2*dynVar(3)*sigma*dynVar(1)/m/dynVar(2)^2+dynVar(4)*e/epsilon, r-sigma*(1-dynVar(1)/m/dynVar(2))+sigma*dynVar(1)/m/dynVar(2)+gamma/(C+tau)*dynVar(1)^theta1*theta1/dynVar(1)/(1+dynVar(1)^theta2)-gamma/(C+tau)*dynVar(1)^theta1/(1+dynVar(1)^theta2)^2*dynVar(1)^theta2*theta2/dynVar(1)+eta*u(1), e*dynVar(2)/epsilon; ...
            -2*dynVar(3)*sigma*dynVar(1)/m/dynVar(2)^2+dynVar(4)*e/epsilon, 2*dynVar(3)*sigma*dynVar(1)^2/m/dynVar(2)^3, -sigma*dynVar(1)^2/m/dynVar(2)^2, r-(-d-e*dynVar(1))/epsilon];
end
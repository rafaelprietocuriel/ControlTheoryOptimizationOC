function out=harvest2DCanonicalSystemIC(t,dynVar,pararg,arcarg)
%
% returns the canonical system plus the algebraic eqaution for the control value of harvest2D model
                                                                                            
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
        u=dynVar(5,:);
        out=[[sigma.*dynVar(1,:).*(1-dynVar(1,:)./m./dynVar(2,:))-gamma./(C+tau).*dynVar(1,:).^theta1./(1+dynVar(1,:).^theta2)-eta.*u(1,:).*dynVar(1,:)]; ...
            [(n-d.*dynVar(2,:)-e.*dynVar(2,:).*dynVar(1,:))./epsilon]; ...
            [r.*dynVar(3,:)-p.*u(1,:)-dynVar(3,:).*(sigma.*(1-dynVar(1,:)./m./dynVar(2,:))-sigma.*dynVar(1,:)./m./dynVar(2,:)-gamma./(C+tau).*dynVar(1,:).^theta1.*theta1./dynVar(1,:)./(1+dynVar(1,:).^theta2)+gamma./(C+tau).*dynVar(1,:).^theta1./(1+dynVar(1,:).^theta2).^2.*dynVar(1,:).^theta2.*theta2./dynVar(1,:)-eta.*u(1,:))+dynVar(4,:).*e.*dynVar(2,:)./epsilon]; ...
            [r.*dynVar(4,:)-dynVar(3,:).*sigma.*dynVar(1,:).^2./m./dynVar(2,:).^2-dynVar(4,:).*(-d-e.*dynVar(1,:))./epsilon]; ...
            p.*dynVar(1,:)-2.*phi.*u(1,:)-dynVar(3,:).*eta.*dynVar(1,:)];
    case 2
        u=harvest2DOptimalControl(t,dynVar,pararg,arcarg);
        out=[[sigma.*dynVar(1,:).*(1-dynVar(1,:)./m./dynVar(2,:))-gamma./(C+tau).*dynVar(1,:).^theta1./(1+dynVar(1,:).^theta2)-eta.*u(1,:).*dynVar(1,:)]; ...
            [(n-d.*dynVar(2,:)-e.*dynVar(2,:).*dynVar(1,:))./epsilon]; ...
            [r.*dynVar(3,:)-p.*u(1,:)-dynVar(3,:).*(sigma.*(1-dynVar(1,:)./m./dynVar(2,:))-sigma.*dynVar(1,:)./m./dynVar(2,:)-gamma./(C+tau).*dynVar(1,:).^theta1.*theta1./dynVar(1,:)./(1+dynVar(1,:).^theta2)+gamma./(C+tau).*dynVar(1,:).^theta1./(1+dynVar(1,:).^theta2).^2.*dynVar(1,:).^theta2.*theta2./dynVar(1,:)-eta.*u(1,:))+dynVar(4,:).*e.*dynVar(2,:)./epsilon]; ...
            [r.*dynVar(4,:)-dynVar(3,:).*sigma.*dynVar(1,:).^2./m./dynVar(2,:).^2-dynVar(4,:).*(-d-e.*dynVar(1,:))./epsilon]];
end

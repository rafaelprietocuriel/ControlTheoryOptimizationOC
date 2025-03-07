function t=absolutetime(ocTrj)

leftarcindex=ocTrj.arcposition(1,:);
rightarcindex=ocTrj.arcposition(2,:);
t=zeros(1,rightarcindex(end));
for arc=1:arcnum(ocTrj)
    diffarctime=diff(ocTrj.arcinterval);
    s=ocTrj.x(leftarcindex(arc):rightarcindex(arc));
    t(leftarcindex(arc):rightarcindex(arc))=diffarctime(arc)*s+(ocTrj.arcinterval(arc)-diffarctime(arc)*(arc-1)); 
end
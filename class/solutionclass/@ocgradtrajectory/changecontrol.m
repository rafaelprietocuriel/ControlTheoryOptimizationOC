function ocgTrj=changecontrol(ocgTrj,ocObj,ocTrj)

t=time(ocObj,ocTrj,1);
v=control(ocObj,ocTrj,1);
[t,idx]=unique(t);
v=v(:,idx);
tnew=time(ocgTrj);
if size(ocgTrj.variable.v,1)==1
    ocgTrj.variable.v=interp1(t,v,tnew);
else
    ocgTrj.variable.v=interp1(t,v.',tnew).';
end
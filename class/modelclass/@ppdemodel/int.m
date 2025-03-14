function out=int(ppdeObj,ppdeTrj,varargin)

t=triangles(ppdeTrj);
p=points(ppdeTrj);

a1=t(1,:);a2=t(2,:);a3=t(3,:);% Corner point indices

% Triangle sides
r23x=p(1,a3)-p(1,a2);r23y=p(2,a3)-p(2,a2);r31x=p(1,a1)-p(1,a3);
r31y=p(2,a1)-p(2,a3);
ar=abs(r31x.*r23y-r31y.*r23x)/2;% areas of tiangles 
ut=pdeintrp(p,t,u);
out=ut*ar';

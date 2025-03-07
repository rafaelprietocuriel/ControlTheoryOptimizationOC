function nrm=norm(ppdeObj,solObj,varargin)
%
% 
nrm=[];
normtype=[];
spec=[];
if isempty(ppdeObj)
    return
end
if nargin<=3
    spec='state';
end
if nargin==2
    normtype=2;
end
if isempty(normtype)
    normtype=2;
end
if isempty(spec)
    spec='state';
end
par=parametervalue(ppdeObj);
if isppdeprimitive(solObj)
    arcarg=0;
    t=[0 1];
    depvar=[solObj.y solObj.y];
    pt=points(solObj);
    eg=edges(solObj);
    tr=triangles(solObj);
    trianglearea=pdetrg(pt,tr);
    epdeSol=solObj;
elseif isppdetrajectory(solObj) || isppdeasymptotic(solObj)
    arcarg=0;
    t=time(ppdeObj,solObj,1);
    depvar=genstates(solObj);
    pt=points(solObj);
    eg=edges(solObj);
    tr=triangles(solObj);
    trianglearea=pdetrg(pt,tr);
    if isppdeasymptotic(solObj)
        epdeSol=ppdeprimitive(solObj);
    else
        epdeSol=ppdeprimitive([]);
    end
end
vol=sqrt((max(pt(1,:))-min(pt(1,:)))*(max(pt(2,:))-min(pt(2,:))));

depvar=depvar(meshindex(ppdeObj,solObj,spec),:);
n=length(t);
tmp=zeros(1,n);
switch normtype
    case 1
        for ii=1:n
            tmp(ii)=sum(trianglearea.*pdeintrp(eg,tr,abs(depvar(:,ii))));
        end
        nrm=sum(diff(t).*(tmp(1:n-1)+tmp(2:n))/2)/vol;
    case 2
        for ii=1:n
            tmp(ii)=sqrt(sum(trianglearea.*pdeintrp(eg,tr,depvar(:,ii)).^2));
        end
        nrm=sum(diff(t).*(tmp(1:n-1)+tmp(2:n))/2)/vol;
    case inf
        nrm=max(depvar(:));
end

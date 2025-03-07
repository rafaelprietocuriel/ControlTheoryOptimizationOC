function nrm=spatialnorm(ppdeObj,solObj,varargin)
%
% 
nrm=[];
spec=[];
if isempty(ppdeObj)
    return
end
if nargin>=3
    spec=varargin{1};
end
if isempty(spec)
    spec='state';
end

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

depvar=depvar(meshindex(ppdeObj,epdeSol,spec),:);

for ii=1:length(t)
    nrm(ii)=sqrt(sum(trianglearea.*pdeintrp(eg,tr,depvar(:,ii)).^2));
end

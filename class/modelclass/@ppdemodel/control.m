function u=control(ppdeObj,solObj,varargin)
%

u=[];
if isempty(ppdeObj) || isempty(solObj)
    return
end
par=parametervalue(ppdeObj);

if isppdeprimitive(solObj) ||  isppdetrajectory(solObj) || isppdeasymptotic(solObj)
    depvar=dependentvar(solObj);
    pt=points(solObj);
    npt=numpoints(solObj);
    t=time(ppdeObj,solObj);
    nt=length(t);
end

coeffidx=meshindex(ppdeObj,solObj,'depvar');
u=zeros(controlnum(ppdeObj)*npt,nt);
for ii=1:nt
    actdepvar=depvar(:,ii);
    tmp=feval(ppdeObj,'OptimalControlDistributed',t(ii),pt,actdepvar(coeffidx),par,0).';
    u(:,ii)=tmp(:);
end


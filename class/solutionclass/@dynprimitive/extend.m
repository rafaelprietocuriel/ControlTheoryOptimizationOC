function extdynPrim=extend(dynPrim,T)
%
extdynPrim=struct(dynPrim);
if isperiodic(dynPrim)
    T=T(end);
    n=ceil(T/period(dynPrim));
    repmat(0:n-1,length(dynPrim.octrajectory.x),1);
    sol.x=repmat(dynPrim.octrajectory.x,1,n)+ans(:)';
    sol.y=repmat(dynPrim.octrajectory.y,1,n);
    removeidx=(1:n-1)*length(dynPrim.octrajectory.x);
    sol.x(removeidx)=[];
    sol.x=sol.x/sol.x(end);
    sol.y(:,removeidx)=[];
    sol.parameters=[];
    sol.arcinterval=[0 n*period(dynPrim)]; %
    sol.timehorizon=inf;
    sol.arcarg=arcargument(dynPrim);
    sol.x0=sol.x(1);
    arcposition=find(diff(sol.x)==0);
    sol.arcposition=[[1 arcposition+1];[arcposition numel(sol.x)]];
    extdynPrim.octrajectory=octrajectory(sol);
    extdynPrim.octrajectory.linearization=linearization(dynPrim);
    extdynPrim=dynprimitive(extdynPrim);
elseif isequilibrium(dynPrim)
    extdynPrim.octrajectory.x=linspace(0,1,length(T)+1); %
    extdynPrim.octrajectory.y=extdynPrim.octrajectory.y(:,ones(1,length(T)+1)); %
    extdynPrim.octrajectory.arcinterval=[0 T(end)]; %
    extdynPrim.octrajectory.arcposition=[1;length(T)+1]; %
    extdynPrim.octrajectory=octrajectory(extdynPrim.octrajectory);
    extdynPrim.octrajectory.linearization=linearization(dynPrim);
    extdynPrim=dynprimitive(extdynPrim);
end

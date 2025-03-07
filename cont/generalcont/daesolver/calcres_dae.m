function res=calcres_dae(tcol,x,y,z1,freepar,modelpar,dae,bc,ic,tangent)
global OCMATCONT

if ~isempty(tangent)
    res=zeros(OCMATCONT.MESHDATA.continuationindex,1);
else
    res=zeros(OCMATCONT.MESHDATA.continuationindex-1,1);
end
res(OCMATCONT.bcidx)=bc(x(:,OCMATCONT.MESHDATA.arcpositioncollocationmesh(1,:)),x(:,OCMATCONT.MESHDATA.arcpositioncollocationmesh(2,:)),freepar,modelpar);

for arc = 1:OCMATCONT.arcnumber
    % retrieves indices of the collocation points (idxmesh) and at the
    % meshpoints (idxbd)
    idxmesh=OCMATCONT.MESHDATA.arcpositioncollocationmesh(1,arc):OCMATCONT.MESHDATA.arcpositioncollocationmesh(2,arc)-1;
    idxbd=idxmesh(1:OCMATCONT.CollocationNumber+1:end);
    idxmesh(idxbd)=[];
    
    % z_ijk-f()=0 at the collocation points the first order components
    % satisfy the ODE
    fz=dae(tcol(idxmesh),x(:,idxmesh),freepar,modelpar);
    res(OCMATCONT.MESHDATA.firstorderFinterioridx(OCMATCONT.MESHDATA.arcpart4firstorderFinterioridx(1,arc):OCMATCONT.MESHDATA.arcpart4firstorderFinterioridx(2,arc)))= ...
        z1(:,OCMATCONT.MESHDATA.arcpart4firstorderzcoefficientidx(1,arc):OCMATCONT.MESHDATA.arcpart4firstorderzcoefficientidx(2,arc))-fz(OCMATCONT.firstordercoordinate,:);
    
    % f()=0 at the collocation points the zero order equations are
    % satisfied
    res(OCMATCONT.MESHDATA.zeroorderFinterioridx(OCMATCONT.MESHDATA.arcpart4zeroorderFinterioridx(1,arc):OCMATCONT.MESHDATA.arcpart4zeroorderFinterioridx(2,arc)))= ...
        fz(OCMATCONT.zeroordercoordinate,:);

    %P(0)-y=0 at the interior mesh points the first order components are
    %continuous
    res(OCMATCONT.MESHDATA.transitionFidx(OCMATCONT.MESHDATA.arcpart4transitionFidx(1,arc):OCMATCONT.MESHDATA.arcpart4transitionFidx(2,arc)))= ...
        x(OCMATCONT.firstordercoordinate,idxbd(2:end))-y(:,(OCMATCONT.MESHDATA.arcpart4ycoefficientidx(1,arc)+1):OCMATCONT.MESHDATA.arcpart4ycoefficientidx(2,arc));
end


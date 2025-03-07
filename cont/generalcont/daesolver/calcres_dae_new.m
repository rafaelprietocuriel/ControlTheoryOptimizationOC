function res=calcres_dae_new(tmesh,x,y,z1,freepar,modelpar,dae,bc,ic,tangent)
global OCMATCONT

N=diff(getarcposition(tmesh));
continuationindex=sum(N)*((OCMATCONT.CollocationNumber+1)*OCMATCONT.firstordercomponentnumber+OCMATCONT.CollocationNumber*OCMATCONT.zeroordercomponentnumber)+OCMATCONT.freeparameternum;
[arcpart4yidx,arcpart4fzidx,arcpart4fFintidx,arcpart4zFintidx,arcpart4tFintidx]=calcresarcpartidx_dae(N);

tcol=makecollocationmesh(tmesh);
arcposcol=getarcpositioncol(tcol);
if ~isempty(tangent)
    res=zeros(continuationindex,1);
else
    res=zeros(continuationindex-1,1);
end
res(OCMATCONT.bcidx)=bc(x(:,arcposcol(1,:)),x(:,arcposcol(2,:)),freepar,modelpar);

for arc = 1:OCMATCONT.arcnumber
    % retrieves indices of the collocation points (idxmesh) and at the
    % meshpoints (idxbd)
    idxmesh=arcposcol(1,arc):arcposcol(2,arc)-1;
    idxbd=idxmesh(1:OCMATCONT.CollocationNumber+1:end);
    idxmesh(idxbd)=[];
    [Iidx,FoIidx,ZoIidx,TIidx]=calcresmatrices_dae(OCMATCONT.order,N(arc));
    
    % z_ijk-f()=0 at the collocation points the first order components
    % satisfy the ODE
    fz=dae(tcol(idxmesh),x(:,idxmesh),freepar,modelpar);
    res(FoIidx(arcpart4fFintidx(1,arc):arcpart4fFintidx(2,arc)))= ...
        z1(:,arcpart4fzidx(1,arc):arcpart4fzidx(2,arc))-fz(OCMATCONT.firstordercoordinate,:);
    
    % f()=0 at the collocation points the zero order equations are
    % satisfied
    res(ZoIidx(arcpart4zFintidx(1,arc):arcpart4zFintidx(2,arc)))= ...
        fz(OCMATCONT.zeroordercoordinate,:);

    %P(0)-y=0 at the interior mesh points the first order components are
    %continuous
    res(TIidx(arcpart4tFintidx(1,arc):arcpart4tFintidx(2,arc)))= ...
        x(OCMATCONT.firstordercoordinate,idxbd(2:end))-y(:,(arcpart4yidx(1,arc)+1):arcpart4yidx(2,arc));
end


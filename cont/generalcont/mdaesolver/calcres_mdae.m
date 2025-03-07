function res=calcres_mdae(t,x,y,z1,freepar,modelpar,dae,bc,ic,tangent)
global OCMATCONT

continuationindex=sum(OCMATCONT.MeshNumber-1)*((OCMATCONT.CollocationNumber+1)*sum(OCMATCONT.firstordercomponentnumber)+OCMATCONT.CollocationNumber*sum(OCMATCONT.zeroordercomponentnumber))+OCMATCONT.freeparameternum;

tcol=makecollocationmesh(t);
arcposcol=getarcpositioncol(tcol);
coeffpos=getcoeffposition(t);
if ~isempty(tangent)
    res=zeros(continuationindex,1);
else
    res=zeros(continuationindex-1,1);
end

ctr=0;
for stage=1:OCMATCONT.stagenumber
    N=OCMATCONT.MeshNumber(ctr+(1:OCMATCONT.arcnumber(stage)));
    ctr=ctr+OCMATCONT.arcnumber(stage);
    res(OCMATCONT.bcidx{stage})=bc(x(:,arcposcol{stage}(1,:)),x(:,arcposcol{stage}(2,:)),freepar{stage},modelpar{stage});
    [Iidx,FoIidx,ZoIidx,TIidx]=calcresmatrices_mdae(OCMATCONT.order{stage},sum(N));
    for arc = 1:OCMATCONT.arcnumber
        % retrieves indices of the collocation points (idxmesh) and at the
        % meshpoints (idxbd)
        idxmesh=arcposcol{stage}(1,arc):arcposcol{stage}(2,arc)-1;
        idxbd=idxmesh(1:OCMATCONT.CollocationNumber+1:end);
        idxmesh(idxbd)=[];
        [arcpart4yidx,arcpart4fzidx,arcpart4fFintidx,arcpart4zFintidx,arcpart4tFintidx]=calcresarcpartidx_dae(N(arc));


        % z_ijk-f()=0 at the collocation points the first order components
        % satisfy the ODE
        fz=dae{stage}(tcol{stage}(idxmesh),x{stage}(:,idxmesh),freepar,modelpar{stage});
        res(coeffpos(FoIidx(arcpart4fFintidx(1,arc):arcpart4fFintidx(2,arc)),stage))= ...
            z1{stage}(:,arcpart4fzidx(1,arc):arcpart4fzidx(2,arc))-fz(OCMATCONT.firstordercoordinate{stage},:);

        % f()=0 at the collocation points the zero order equations are
        % satisfied
        res(coeffpos(ZoIidx(arcpart4zFintidx(1,arc):arcpart4zFintidx(2,arc)),stage))= ...
            fz(OCMATCONT.zeroordercoordinate{stage},:);

        %P(0)-y=0 at the interior mesh points the first order components are
        %continuous
        res(coeffpos(TIidx(arcpart4tFintidx(1,arc):arcpart4tFintidx(2,arc)),stage))= ...
            x{stage}(OCMATCONT.firstordercoordinate,idxbd(2:end))-y{stage}(:,(arcpart4yidx(1,arc)+1):arcpart4yidx(2,arc));
    end
end

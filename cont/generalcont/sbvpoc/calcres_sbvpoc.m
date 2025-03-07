function res=calcres_sbvpoc(tmesh,x,y,z1,freepar,modelpar,dae,bc)
global OCMATCONT OCBVP

arcnumber=OCBVP.arcnumber;
arcpositioncollocationmesh=OCBVP.arcpositioncollocationmesh;

res=zeros(OCMATCONT.continuationindex,1);
%res(1:nBCs) = bc(tmesh,y(:,1,Lidx),y(:,1,Ridx),z(:,:,Lidx),z(:,:,Ridx),freepar,modelpar);

for arc = 1:arcnumber
    idx=arcpositioncollocationmesh(1,arc):arcpositioncollocationmesh(2,arc);
    fz=dae(tmesh(idx),x(:,OCBVP.firstorderidx,:,idx),freepar,modelpar,arc);
    z_arc=z1(:,idx);
    %Evaluation of the equations for the solver
    fz=zeros(numeq,numcols,2);
    for jj=1:nummeshintv(region)
        % dynamics and algebraic equations have to be satisfied at the
        % collocation points
        fz(odecoord,numcolscoord,1)=evaluatePcols(jj,diffmesh_arc,y_arc,z_arc,odecoord,psival,numcols);
        fz(odecoord,numcolscoord,2)=z_arc(odecoord,numcolscoord,jj);
        fz(aecoord,numcolscoord,1)=z_arc(aecoord,numcolscoord,jj);
        fz(:,:,2)-ode(cols(jj,:),fz(:,:,1),FcnArgs{:});
        res(countercols+domainddata(arcindex).numeqnumcolscoord,1)=ans(:);
        countercols=countercols+domainddata(arcindex).numeqnumcols;
        if jj<nummeshintv(region)
            % at the time mesh the ODE solution path has to be continuous
            res(countercols+odecoord)=evaluatePleft(jj,diffmesh_arc,y_arc,z_arc,odecoord,psival,numcols)-y_arc(odecoord,1,jj+1);
            countercols=countercols+domainddata(arcindex).numode;
        end
    end
end
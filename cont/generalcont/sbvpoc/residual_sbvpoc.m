function [res resindex]=residual_sbvpoc(tmesh,y,z,freepar,contval,modelpar,odefile,respoints)
global OCMATCONT
domainddata=OCMATCONT.DOMAINDDATA;

leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
diffmesh=diff(tmesh);
counter=0;
resindex=ones(OCMATCONT.HE.numarc,2);
%Determine the interval where the residualpoint is located
for arc=1:OCMATCONT.HE.numarc
    arcidx=leftarcindex(arc)-arc+1:rightarcindex(arc)-arc;
    arcindex=OCMATCONT.HE.arcindex(arc); % the arcindex of the actual arc
    eqcoord=domainddata(arcindex).eqcoord;
    odecoord=domainddata(arcindex).odecoord;
    aecoord=domainddata(arcindex).aecoord;
    numeq=domainddata(arcindex).numeq;
    numode=domainddata(arcindex).numode;
    numae=domainddata(arcindex).numae;
    numcols=domainddata(arcindex).numcols;
    numcolscoord=domainddata(arcindex).numcolscoord;
    psi=domainddata(arcindex).psi;
    psi0=domainddata(arcindex).psi0;
    y_arc=y(odecoord,:,arcidx);
    z_arc=z(eqcoord,numcolscoord,arcidx);
    diffmesh_arc=diffmesh(leftarcindex(arc):rightarcindex(arc)-1);
    tmesh_arc=tmesh(leftarcindex(arc):rightarcindex(arc));
    %Evaluation of the equations for the solver
    fz=zeros(numeq,2);
    while ~isempty(respoints) && respoints(1)<tmesh(rightarcindex(arc))
        counter=counter+1;
        %Determine values of polynomials of the residualpoint in the
        %corresponding interval
        idx=find(respoints(1)<tmesh_arc);
        idx=idx(1)-1;
        psival=[];
        psival0=[];
        for jj=numcolscoord
            psival(1,jj,1)=polyval(psi(jj,:,1),(respoints(1)-tmesh_arc(idx))/diffmesh_arc(idx));
            psival0(1,jj,1)=polyval(psi0(jj,:,1),(respoints(1)-tmesh_arc(idx))/diffmesh_arc(idx));
        end
        fz(odecoord,1)=y_arc(odecoord,1,idx)+sum(diffmesh_arc(idx)*z_arc(odecoord,:,idx).*psival(ones(numode,1),:,1),2);
        fz(odecoord,2)=sum(z_arc(odecoord,:,idx).*psival0(ones(numode,1),:,1),2);
        if numae
            fz(aecoord,1)=sum(z_arc(aecoord,:,idx).*psival0(ones(numae,1),:,1),2);
        end
        res(eqcoord,counter)=fz(:,2)-odefile(respoints(1),fz(:,1),modelpar,freepar,contval,arc);
        respoints(1)=[];
    end
    resindex(arc,2)=counter;
    if arc<OCMATCONT.HE.numarc
        resindex(arc+1,1)=counter+1;
    end
end

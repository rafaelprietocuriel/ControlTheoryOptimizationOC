function sol=evalatmesh(tmesh,y,z)

global OCMATCONT OCBVP

switch OCMATCONT.bvpmethod
    case 'sbvpoc'
        domainddata=OCMATCONT.DOMAINDDATA;
        diffmesh=diff(tmesh);

        leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
        rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
        nummeshintv=OCMATCONT.HE.TIMEDDATA.nummeshintv;
        sol.x=tmesh;
        sol.y=zeros(max([domainddata(OCMATCONT.HE.arcindex).numeq]),rightarcindex(OCMATCONT.HE.numarc));
        for arc=1:OCMATCONT.HE.numarc
            arcindex=OCMATCONT.HE.arcindex(arc);
            numcols=domainddata(arcindex).numcols;
            numcolscoord=domainddata(arcindex).numcolscoord;
            numode=domainddata(arcindex).numode;
            numae=domainddata(arcindex).numae;
            psival=domainddata(arcindex).psival;
            odecoord=domainddata(arcindex).odecoord; % coordinate of (first order) odes
            diffmesh_arc=diffmesh(leftarcindex(arc):rightarcindex(arc)-1);
            y_arc=y(odecoord,:,leftarcindex(arc)-arc+1:rightarcindex(arc)-arc);
            z_arc=z(:,numcolscoord,leftarcindex(arc)-arc+1:rightarcindex(arc)-arc);
            sol.y(odecoord,leftarcindex(arc):rightarcindex(arc)-1)=y_arc;
            sol.y(odecoord,rightarcindex(arc))=y_arc(odecoord,:,nummeshintv(arc))+sum(diffmesh_arc(nummeshintv(arc))*z_arc(odecoord,:,nummeshintv(arc)).*psival(ones(1,numode),:,numcols+1),2);
            if numae
                aecoord=domainddata(arcindex).aecoord; % coordinate of algebraic equations
                psival0=domainddata(arcindex).psival0;
                pos=numcols+2;
                sol.y(aecoord,leftarcindex(arc):rightarcindex(arc)-1)=sum(z_arc(aecoord,:,:).*psival0(ones(1,numae),:,pos(ones(1,nummeshintv(arc)))),2);
                sol.y(aecoord,rightarcindex(arc))=sum(z_arc(aecoord,:,nummeshintv(arc)).*psival0(ones(1,numae),:,numcols+1),2);
            end
        end
    case {'bvp4c','gbvp4c','bvp6c','bvp5c'}
        sol.x=tmesh;
        sol.y=y;
end

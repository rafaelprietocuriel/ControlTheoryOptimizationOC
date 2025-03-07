function [yl,yr,diffmeshatswitch]=boundarycollocationpoints(tmesh,ya,yb,za,zb)
global OCMATCONT OCBVP
domainddata=OCMATCONT.DOMAINDDATA;
diffmesh=diff(tmesh);
yl=zeros(OCMATCONT.HE.maxnumeq,OCMATCONT.HE.numarc);
yr=zeros(OCMATCONT.HE.maxnumeq,OCMATCONT.HE.numarc);
diffmeshatswitch=zeros(OCMATCONT.HE.numarc,4); % determines the mesh size adjacent to the arc time intervals

leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;
for arc=1:OCMATCONT.HE.numarc
    % we have to distinguish between the arc index of the time mesh
    % and the values of the coefficients, for the natvie MATLAB
    % solvers these two indices coincide. For SBVPOC the right side
    % values are only given by the interpolation formula
    tmeshidx=[leftarcindex(arc) rightarcindex(arc)-1];
    idx=[OCBVP.Lidx(arc) OCBVP.Ridx(arc)];
    diffmeshatswitch(arc,2:3)=diffmesh(tmeshidx);
    if arc>1
        diffmeshatswitch(arc,1)=diffmesh(tmeshidx(1)-2);
    end
    if arc<OCMATCONT.HE.numarc
        diffmeshatswitch(arc,4)=diffmesh(tmeshidx(2)+2);
    end
    arcindex=OCMATCONT.HE.arcindex(arc);
    aecoord=domainddata(arcindex).aecoord;
    odecoord=domainddata(arcindex).odecoord;
    numae=domainddata(arcindex).numae;
    numode=domainddata(arcindex).numode;
    numcolscoord=domainddata(arcindex).numcolscoord;
    numcols=domainddata(arcindex).numcols;

    yl(odecoord,arc)=ya(odecoord,1,arc);
 %               yl(domainddata(arcindex).odecoord,arc)=evalcollpoly(idx(1),diffmesh,y,z,arc,0,0);
%                 if abs(yl(domainddata(arcindex).odecoord,arc)-y(odecoord,1,idx(1)))>1e-12
%                     pause
%                 end
    if ~isempty(aecoord)
        %yl(aecoord,arc)=sum(z(aecoord,numcolscoord,idx(1)).*domainddata(arcindex).psival0(ones(numae,1),:,numcols+2),2);
                        yl(aecoord,arc)=evalcollpoly(arc,diffmesh,ya,za,arc,0,1);
%                         if abs(yl(domainddata(arcindex).aecoord,arc)-(sum(z(aecoord,numcolscoord,idx(1)).*domainddata(arcindex).psival0(ones(numae,1),:,numcols+2),2)))>1e-12
%                             pause
%                         end
    end
    yr(odecoord,arc)=yb(odecoord,1,arc)+sum(diffmesh(rightarcindex(arc)-1)*zb(odecoord,numcolscoord,arc).*domainddata(arcindex).psival(ones(numode,1),:,numcols+1),2);
                %yr(odecoord,arc)=evalcollpoly(arc,diffmesh,yb,zb,arc,1,0);
%                 if abs(yr(domainddata(arcindex).odecoord,arc)-(y(odecoord,1,idx(2))+sum(diffmesh(idx(2))*z(odecoord,numcolscoord,idx(2)).*domainddata(arcindex).psival(ones(numode,1),:,numcols+1),2)))>1e-12
%                     pause
%                 end
    if ~isempty(domainddata(arcindex).aecoord)
        yr(aecoord,arc)=sum(zb(aecoord,numcolscoord,arc).*domainddata(arcindex).psival0(ones(numae,1),:,numcols+1),2);
        yr(domainddata(arcindex).aecoord,arc)=evalcollpoly(arc,diffmesh,yb,zb,arc,1,1);
%                         if abs(yr(domainddata(arcindex).aecoord,arc)-(sum(z(aecoord,numcolscoord,idx(2)).*domainddata(arcindex).psival0(ones(numae,1),:,numcols+1),2)))>1e-12
%                             pause
%                         end
    end
end

OCMATCONT.HE.yL=yl; % L means left side interval
OCMATCONT.HE.yR=yr; % R means right side interval
OCMATCONT.HE.gridatswitch=diffmeshatswitch;
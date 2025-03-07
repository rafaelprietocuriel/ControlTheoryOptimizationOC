function [yal yar ybl ybr diffmeshatswitch]=determineswitchstate(tmesh,y,z,arc)
global OCMATCONT
ya=OCMATCONT.HE.yL; % L means left side interval
yb=OCMATCONT.HE.yR; % R means right side interval
diffmeshatswitch=[];

if arc<OCMATCONT.HE.numarc
    if arc==1
        yal=[];
        yar=ya(:,arc);
        ybl=yb(:,arc);
        ybr=ya(:,arc+1);
    else
        yal=ya(:,arc-1);
        yar=yb(:,arc);
        ybl=yb(:,arc);
        ybr=yb(:,arc+1);
    end
end
if OCMATCONT.HE.numarc==1
    yal=[];
    yar=ya(:,arc);
    ybl=yb(:,arc);
    ybr=[];
elseif arc==OCMATCONT.HE.numarc
    yal=ya(:,arc-1);
    yar=ya(:,arc);
    ybl=yb(:,arc);
    ybr=[];
end

if strcmp(OCMATCONT.solver,'sbvpoc')
    diffmeshatswitch=OCMATCONT.HE.gridatswitch(arc,1:4);
end
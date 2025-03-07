function [tcol,tmeshidx]=makecollocationmesh(tmesh)
global OCMATCONT

% returns the mesh + collocation points

arcposition=getarcposition(tmesh);
tcol=zeros(1,(length(tmesh)-OCMATCONT.arcnumber)*(OCMATCONT.CollocationNumber+1));
tmeshidx=zeros(1,OCMATCONT.arcnumber*length(tmesh));
ctr=0;
ctr1=0;
for ii=1:OCMATCONT.arcnumber
    tarc=tmesh(arcposition(1,ii):arcposition(2,ii));
    Narc=length(tarc);
    h=diff(tarc);
    tmptfine=[tarc(1:Narc-1);myrepmat(tarc(1:Narc-1),OCMATCONT.CollocationNumber,1,1,Narc-1)+myrepmat(h,OCMATCONT.CollocationNumber,1,1,Narc-1).*myrepmat(OCMATCONT.CollocationPoint(:),1,Narc-1,OCMATCONT.CollocationNumber,1)];
    tcol(ctr+(1:(Narc-1)*(OCMATCONT.CollocationNumber+1)+1))=[tmptfine(:).' tarc(Narc)];
    ctr=ctr+(Narc-1)*(OCMATCONT.CollocationNumber+1)+1;
    tmeshidx(ctr1+(1:Narc))=ctr1+(1:OCMATCONT.CollocationNumber+1:Narc*(OCMATCONT.CollocationNumber+1));
    ctr1=ctr1+Narc;
end

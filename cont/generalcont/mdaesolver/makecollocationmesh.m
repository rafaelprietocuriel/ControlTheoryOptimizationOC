function tcol=makecollocationmesh(tmesh)
global OCMATCONT

% returns the mesh + collocation points

arcposition=getarcposition(tmesh);
tcol=cell(1,OCMATCONT.stagenumber);
for ii=1:OCMATCONT.stagenumber
    tcol{ii}=zeros(1,(length(tmesh{ii})-OCMATCONT.arcnumber(ii))*(OCMATCONT.CollocationNumber+1));
    ctr=0;
    for jj=1:OCMATCONT.arcnumber(ii)
        tarc=tmesh(arcposition{ii}(1,jj):arcposition{ii}(2,jj));
        Narc=length(tarc);
        h=diff(tarc);
        tmptfine=[tarc(1:Narc-1);myrepmat(tarc(1:Narc-1),OCMATCONT.CollocationNumber,1,1,Narc-1)+myrepmat(h,OCMATCONT.CollocationNumber,1,1,Narc-1).*myrepmat(OCMATCONT.CollocationPoint(:),1,Narc-1,OCMATCONT.CollocationNumber,1)];
        tcol{ii}(ctr+(1:(Narc-1)*(OCMATCONT.CollocationNumber+1)+1))=[tmptfine(:).' tarc(Narc)];
        ctr=ctr+(Narc-1)*(OCMATCONT.CollocationNumber+1)+1;
    end
end

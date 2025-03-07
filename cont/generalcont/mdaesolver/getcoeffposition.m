function coeffpos=getcoeffposition(t)
global OCMATCONT
% returns the position of the arcs for the collocation mesh for each stage

arcpos=getarcposition(t);
ctr=OCMATCONT.bcnum;
coeffposition=zeos(2,OCMATCONT.stagenumber);
for ii=1:OCMATCONT.stagenumber
    coeffposition(:,ii)=[ctr+1;arcpos{ii}(2,end) (arcpos{ii}(2,end)-1)*((OCMATCONT.CollocationNumber+1)*OCMATCONT.firstordercomponentnumber(ii)+OCMATCONT.CollocationNumber*OCMATCONT.zeroordercomponentnumber(ii))];
end
coeffpos=coeffposition;

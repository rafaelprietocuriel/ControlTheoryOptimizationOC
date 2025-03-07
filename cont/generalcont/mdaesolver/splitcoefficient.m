function [coeff1,coeff0]=splitcoefficient(coeff0,Narc,arc)
% the input argument coeff0 is the total coefficient, the output argument
% coeff0, are the coefficients of the zero order components

global OCMATCONT
persistent foi zoi

if isempty(foi) || length(foi)~=OCMATCONT.zeroordercomponentnumber*(Narc-1)
    foi=myrepmat(OCMATCONT.firstordercoordinate.',1,Narc-1,1,OCMATCONT.firstordercomponentnumber)+reshape(myrepmat((0:Narc-2)*OCMATCONT.componentnumber,OCMATCONT.firstordercomponentnumber,1,1,Narc-1),1,[]);
    zoi=myrepmat(OCMATCONT.zeroordercoordinate.',1,Narc-1,1,OCMATCONT.zeroordercomponentnumber)+reshape(myrepmat((0:Narc-2)*OCMATCONT.componentnumber,OCMATCONT.zeroordercomponentnumber,1,1,Narc-1),1,[]);
end
coeffidx=myrepmat((1:OCMATCONT.componentnumber*OCMATCONT.CollocationNumber+OCMATCONT.firstordercomponentnumber).',1,Narc-1,OCMATCONT.componentnumber*OCMATCONT.CollocationNumber+OCMATCONT.firstordercomponentnumber,1)+ ...
    myrepmat((0:Narc-2)*(OCMATCONT.firstordercomponentnumber+OCMATCONT.componentnumber*OCMATCONT.CollocationNumber),OCMATCONT.componentnumber*OCMATCONT.CollocationNumber+OCMATCONT.firstordercomponentnumber,1,1,Narc-1);
ycoefficientidx=coeffidx(1:OCMATCONT.firstordercomponentnumber,:);
Z=coeffidx(OCMATCONT.firstordercomponentnumber+(1:OCMATCONT.CollocationNumber*OCMATCONT.componentnumber),:);
Z=reshape(Z,OCMATCONT.componentnumber,[]);
zcoefficientidx=zeros((Narc-1)*OCMATCONT.componentnumber,OCMATCONT.CollocationNumber);
tmp=zeros(OCMATCONT.componentnumber,OCMATCONT.CollocationNumber);
for ii=1:Narc-1
    tmp(:)=Z(:,(ii-1)*OCMATCONT.CollocationNumber+1:ii*OCMATCONT.CollocationNumber);
    zcoefficientidx((ii-1)*OCMATCONT.componentnumber+1:ii*OCMATCONT.componentnumber,:)=tmp;
end
coeff1=coeff0([ycoefficientidx(:).';zcoefficientidx(foi,:).']);
coeff0=coeff0(zcoefficientidx(zoi,:).');

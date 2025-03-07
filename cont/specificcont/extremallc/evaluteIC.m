function out=evaluteIC(y,z,zref,ocMFVar,limitcycleData,N,gridsize)

% evaluate integral constraints for phase condition

h2=gridsize.*gridsize;
z=permute(z,[1 3 2]);
zref=permute(zref,[1 3 2]);
y=permute(y,[1 3 2]);
zhat=limitcycleData.template.zerosnN;
ztilde=limitcycleData.template.zerosnNncol;
for ii=ocMFVar.ncolcoord
    zreftw=zref(:,:,ii)*limitcycleData.gw(ii);
    zhat=zhat+zreftw;
    zreftwm(:,:,ii)=zreftw;
end
for ii=ocMFVar.ncolcoord
    ztilde(:,:,ii)=sum(limitcycleData.psival3(limitcycleData.template.onesn,ii+limitcycleData.template.zerosN,:).*z(:,:,ii+limitcycleData.template.zerosncol).*zreftwm,3);
end
out=gridsize(limitcycleData.template.onesn,:).*zhat.*y+h2(limitcycleData.template.onesn,:).*sum(ztilde,3);
out=sum(out(:));

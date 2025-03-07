function DI=evaluteDIC(y,z,zref,limitcycleData,ocMFVar,N,gridsize)

% evaluate integral constraints for phase condition

h2=gridsize.*gridsize;
zref=permute(zref,[1 3 2]);
DIy_=limitcycleData.template.zerosnN;
DIz_=limitcycleData.template.zerosnNncol;
zreftwm=limitcycleData.template.zerosnNncol;
for ii=ocMFVar.ncolcoord
    zreftw=zref(:,:,ii)*limitcycleData.gw(ii);
    DIy_=DIy_+zreftw;
    zreftwm(:,:,ii)=zreftw;
end
for ii=ocMFVar.ncolcoord
        DIz_(:,:,ii)=sum(limitcycleData.psival3(limitcycleData.template.onesn,ii+limitcycleData.template.zerosN,:).*zreftwm,3);
end
DIy_=gridsize(limitcycleData.template.onesn,:).*DIy_;
DIz_=h2(limitcycleData.template.onesn,:,limitcycleData.template.onesncol).*DIz_;
DI=cat(3,DIy_,DIz_);
DI=permute(DI,[1 3 2]);
%DI=DI(:,:).';
DI=[DI(:).' limitcycleData.template.zeros1numpar];

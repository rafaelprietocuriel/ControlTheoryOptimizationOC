function [Iidx,FoIidx,ZoIidx,TIidx]=calcresmatrices_dae(order,Nm1)
% 
% for the computation the solution components are reordered such that the
% first order components are at the beginnig and thereafter the zero order
% components. The sol.y components are in the initial ordering.

global OCMATCONT
persistent N germ csgerm

if isempty(N) || N~=Nm1
    N=Nm1;
    germ=order(:,ones(OCMATCONT.CollocationNumber,1));
    germ=[germ(:);OCMATCONT.firstordercomponentnumber];
    csgerm=cumsum([ones(OCMATCONT.componentnumber*OCMATCONT.CollocationNumber,1);OCMATCONT.firstordercomponentnumber]);
    csgerm=csgerm(:,ones(1,N));
end
Iidx=csgerm+myrepmat(cumsum(csgerm(end,:))-csgerm(end,1)+OCMATCONT.bcidx(end),OCMATCONT.componentnumber*OCMATCONT.CollocationNumber+1,1,1,N);

FoIidx=Iidx;FoIidx(germ==0,:)=[];FoIidx(end,:)=[];
FoIidx=FoIidx(:);
ZoIidx=Iidx;ZoIidx(germ==1,:)=[];ZoIidx(end,:)=[];
ZoIidx=ZoIidx(:);

% equations at the mesh points
TIidx=OCMATCONT.bcnum+OCMATCONT.componentnumber*OCMATCONT.CollocationNumber+cumsum(repmat([ones(1,OCMATCONT.firstordercomponentnumber) OCMATCONT.componentnumber*OCMATCONT.CollocationNumber],1,N-1)).';
TIidx(OCMATCONT.firstordercomponentnumber+1:OCMATCONT.firstordercomponentnumber+1:(OCMATCONT.firstordercomponentnumber+1)*(N-1))=[];

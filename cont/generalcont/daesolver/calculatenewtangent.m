function tangent=calculatenewtangent(tmesh,coeff,direction)
global OCMATCONT
EN=[];
EN(OCMATCONT.MESHDATA.continuationindex,1)=1;

DF=OCMATCONT.frechetder(tmesh,coeff,EN,OCMATCONT.dae,OCMATCONT.bc,OCMATCONT.icfun,OCMATCONT.daejac,OCMATCONT.bcjac,OCMATCONT.icfunjac);
%DF(OCMATCONT.MESHDATA.continuationindex,:)=tangent(:).';
[L,U]=lu(DF);            % LU-decomposition of DF


% find the Newton direction
tangent=U\(L\EN);
tangent=sign(tangent(OCMATCONT.MESHDATA.continuationindex))*tangent*direction;
tangent=tangent/norm(tangent);
function [tmesh,coeff,tangent]=adapt2solver(solinit,ode,odejac,bc,bcjac,icfun,icfunjac,odehess,bchess,odetensor3,bctensor3)

clear global OCBVP
global OCMATCONT OCBVP


clear global OCBVP
global OCMATCONT 


OCMATCONT.dae=daefinal;
OCMATCONT.bc=bcfinal;
OCMATCONT.icfun=icfunfinal;
OCMATCONT.daejac=jacfinal;
OCMATCONT.bcjac=bcjacfinal;
OCMATCONT.icfunjac=icfunjacfinal;
OCMATCONT.daehess=daehessfinal;
OCMATCONT.bchess=bchessfinal;
OCMATCONT.daetensor3=daetensor3final;
OCMATCONT.bctensor3=bctensor3final;

        tmesh=solinit.x;
        coeff=solinit.y(:);

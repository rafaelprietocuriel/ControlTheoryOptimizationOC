function [dBCdya,dBCdyb,nCalls,dBCdpar]=gBCnumjac(bc,ya,yb,ExtraArgs)
%BCNUMJAC  Numerically compute dBC/dya, dBC/dyb, and dBC/dpar, if needed.
global OCBVP
dBCoptions=OCBVP.dBCoptions;

if OCBVP.multipointbvp
    % make ya and yb columns vectors
    bcArgs={ya(:),yb(:),OCBVP.maxnumode,bc,ExtraArgs{:}};
    dBCoptions.diffvar=1;
    [dBCdya,ignored,nbc]=numjaccsd(@bcaux,bcArgs,OCBVP.nBCs,dBCoptions);
    nCalls=nbc;
    dBCoptions.diffvar=2;
    [dBCdyb,ignored,nbc]=numjaccsd(@bcaux,bcArgs,OCBVP.nBCs,dBCoptions);
    nCalls=nCalls + nbc;
else
    bcArgs={ya,yb,ExtraArgs{:}};
    dBCoptions.diffvar=1;
    [dBCdya,ignored,nbc]=numjaccsd(bc,bcArgs,OCBVP.nBCs,dBCoptions);
    nCalls=nbc;
    dBCoptions.diffvar=2;
    [dBCdyb,ignored,nbc]=numjaccsd(bc,bcArgs,OCBVP.nBCs,dBCoptions);
    nCalls=nCalls + nbc;
end
if OCBVP.explicitparameterdependence
    bcArgs={ya,yb,ExtraArgs{:}};
    dBCoptions.diffvar=3;
    [dBCdpar,ignored,nbc]=numjaccsd(bc,bcArgs,OCBVP.nBCs,dBCoptions);
    nCalls=nCalls + nbc;
else
    dBCdpar=[];
end
dBCdya=dBCdya(:,OCBVP.keepbccolumns);
dBCdyb=dBCdyb(:,OCBVP.keepbccolumns);

%---------------------------------------------------------------------------

function res=bcaux(Ya,Yb,n,bcfun,varargin)
ya=reshape(Ya,n,[]);
yb=reshape(Yb,n,[]);
res=bcfun(ya,yb,varargin{:});
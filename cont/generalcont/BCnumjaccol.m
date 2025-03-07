function [dBCdya,dBCdyb,nCalls,dBCdpar]=BCnumjaccol(bc,tmesh,ya,yb,za,zb,ExtraArgs)
%BCNUMJAC  Numerically compute dBC/dya, dBC/dyb, and dBC/dpar, if needed.
global OCBVP
dBCoptions=OCBVP.dBCoptions;

% make ya and yb columns vectors
bcArgs={tmesh,ya,yb,za,zb,ExtraArgs{:}};
dBCoptions.diffvar=1;
[dBCdya,ignored,nbc]=numjaccsd(bc,bcArgs,OCBVP.nBCs,dBCoptions);
nCalls=nbc;
dBCoptions.diffvar=2;
[dBCdyb,ignored,nbc]=numjaccsd(bc,bcArgs,OCBVP.nBCs,dBCoptions);
nCalls=nCalls + nbc;
if OCBVP.explicitparameterdependence
    dBCoptions.diffvar=3;
    [dBCdpar,ignored,nbc]=numjaccsd(bc,bcArgs,OCBVP.nBCs,dBCoptions);
    nCalls=nCalls + nbc;
else
    dBCdpar=[];
end

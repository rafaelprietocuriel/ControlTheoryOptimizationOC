function [dBCdya,dBCdyb,nCalls,dBCdpar]=BCnumjac(bc,ya,yb,ExtraArgs)
%BCNUMJAC  Numerically compute dBC/dya, dBC/dyb, and dBC/dpar, if needed.
global OCBVP
dBCoptions=OCBVP.dBCoptions;

if OCBVP.multipointbvp
    % make ya and yb columns vectors
%     bcArgs={ya(:),yb(:),OCBVP.neqn,bc,ExtraArgs{:}};
%     dBCoptions.diffvar=1;
%     [dBCdya,ignored,nbc]=numjaccsd(@bcaux,bcArgs,OCBVP.nBCs,dBCoptions);
%     nCalls=nbc;
%     dBCoptions.diffvar=2;
%     [dBCdyb,ignored,nbc]=numjaccsd(@bcaux,bcArgs,OCBVP.nBCs,dBCoptions);
%     nCalls=nCalls + nbc;


    % Test
    bcArgs = {ya(:),yb(:),OCBVP.neqn,bc,ExtraArgs{:}};
    dBCoptions.thresh = repmat(1e-6,length(ya(:)),1);
    dBCoptions.fac = [];
    dBCoptions.vectvars = []; % BC functions not vectorized

    bcVal  = bcaux(bcArgs{:});
    nCalls = 1;

    dBCoptions.diffvar = 1;
    [dBCdya,ignored,ignored1,nbc] = odenumjac(@bcaux,bcArgs,bcVal,dBCoptions);
    nCalls = nCalls + nbc;
    dBCoptions.diffvar = 2;
    [dBCdyb,ignored,ignored1,nbc] = odenumjac(@bcaux,bcArgs,bcVal,dBCoptions);
    nCalls = nCalls + nbc;
else
%     bcArgs={ya,yb,ExtraArgs{:}};
%     dBCoptions.diffvar=1;
%     [dBCdya,ignored,nbc]=numjaccsd(bc,bcArgs,OCBVP.nBCs,dBCoptions);
%     nCalls=nbc;
%     dBCoptions.diffvar=2;
%     [dBCdyb,ignored,nbc]=numjaccsd(bc,bcArgs,OCBVP.nBCs,dBCoptions);
%     nCalls=nCalls + nbc;

    % Test
    bcArgs = {ya,yb,OCBVP.neqn,bc,ExtraArgs{:}};
    dBCoptions.thresh = repmat(1e-6,length(ya),1);
    dBCoptions.fac = [];
    dBCoptions.vectvars = []; % BC functions not vectorized

    bcVal  = bcaux(bcArgs{:});
    nCalls = 1;

    dBCoptions.diffvar = 1;
    [dBCdya,ignored,ignored1,nbc] = odenumjac(@bcaux,bcArgs,bcVal,dBCoptions);
    nCalls = nCalls + nbc;
    dBCoptions.diffvar = 2;
    [dBCdyb,ignored,ignored1,nbc] = odenumjac(@bcaux,bcArgs,bcVal,dBCoptions);
    nCalls = nCalls + nbc;

end
if OCBVP.explicitparameterdependence
%     bcArgs={ya,yb,ExtraArgs{:}};
%     dBCoptions.diffvar=3;
%     [dBCdpar,ignored,nbc]=numjaccsd(bc,bcArgs,OCBVP.nBCs,dBCoptions);
%     nCalls=nCalls + nbc;

    % Test
    bcArgs = {ya,yb,ExtraArgs{:}};
    dBCoptions.thresh = repmat(1e-6,OCBVP.npar,1);
    dBCoptions.diffvar = 3;
    [dBCdpar,ignored,ignored1,nbc] = odenumjac(bc,bcArgs,bcVal,dBCoptions);
    nCalls = nCalls + nbc;

else
    dBCdpar=[];
end


%---------------------------------------------------------------------------

function res=bcaux(Ya,Yb,n,bcfun,varargin)
ya=reshape(Ya,n,[]);
yb=reshape(Yb,n,[]);
res=bcfun(ya,yb,varargin{:});
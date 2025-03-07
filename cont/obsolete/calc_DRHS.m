function res=calc_DRHS(x,y,freepar,modelpar,map,bc,ic)

global OCMATCONT OCBVP

FcnArgs={0,freepar,modelpar};   

F=zeros(OCBVP.nummap,OCBVP.N-1);
res=zeros(OCMATCONT.HE.numdvariablesmc,1); 

% Boundary conditions
res(1:OCBVP.nBCs)=bc(y(:,OCBVP.Lidx),y(:,OCBVP.Ridx),freepar,modelpar);
phiptr=OCBVP.nBCs;    % active region of res

for arc=1:OCBVP.numarc
    FcnArgs{1}=arc;
    xidx=OCBVP.Lidx(arc):OCBVP.Ridx(arc);   % mesh point index
    xreg=x(xidx(2:end));
    yreg=y(:,xidx);
    Nreg=OCBVP.Nint(arc);

    % differences
    if OCMATCONT.OPTIONS.xyvectorized
        Freg=map(xreg,yreg,FcnArgs{:});
    else
        Freg=zeros(OCBVP.nummap,Nreg);
        for ii=1:Nreg
            Freg(:,ii)=map(xreg(ii),yreg(:,[ii ii+1]),FcnArgs{:});
        end
    end
    res(phiptr+1:phiptr+OCBVP.nummap*OCBVP.Nint(arc))=Freg(:);
    phiptr=phiptr + OCBVP.nummap*OCBVP.Nint(arc);
    F(:,xidx(1:end-1))=Freg;

end
if OCBVP.sumconstraint
    res(phiptr+1:phiptr+OCBVP.sumconstraint)=ic(x,y,freepar,modelpar);
    OCBVP.IC=res(phiptr+1:phiptr+OCBVP.sumconstraint);
end
OCBVP.F=F;



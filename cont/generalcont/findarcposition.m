function findarcposition(tmesh,coeff,tangent,varargin)

global OCMATCONT OCBVP

if length(OCMATCONT.HE.arcarg)==1
    return
end

[failed,infoS]=OCMATCONT.testadmissibility(tmesh,coeff,tangent);

if failed
    infoS
end

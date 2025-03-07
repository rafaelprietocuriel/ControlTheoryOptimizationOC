function failed=DataAdaptation(tmesh,coeff,tangent,tmeshold)
global OCMATCONT

failed=0;
switch OCMATCONT.bvpmethod
    case {'bvp4c'}
        dataadaptation_bvp4c(tmesh);
    case {'bvp5c'}
        dataadaptation_bvp5c(tmesh);
    case {'bvp6c'}
        dataadaptation_bvp6c(tmesh);
    case {'gbvp4c'}
        dataadaptation_gbvp4c(tmesh);
end

if ~isempty(OCMATCONT.dataadaptation) && nargin>1
    failed=OCMATCONT.dataadaptation(tmesh,coeff,tangent,tmeshold);
end
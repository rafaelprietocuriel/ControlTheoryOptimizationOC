function out=arcargument(ocgLim)

out=[];
solInfo=solverinfo(ocgLim);

if isfield(solInfo,'arcarg')
    out=solInfo.arcarg;
end

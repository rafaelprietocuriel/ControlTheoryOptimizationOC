function v=tangent(ocMultPath,varargin)

if isfield(ocMultPath.solverinfo,'tangent')
    v=ocMultPath.solverinfo.tangent;
else
    v=[];
end
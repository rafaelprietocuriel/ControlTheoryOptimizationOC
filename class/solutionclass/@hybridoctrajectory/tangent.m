function v=tangent(ocTrj,varargin)


if isfield(ocTrj.solverinfo,'tangent')
    v=ocTrj.solverinfo.tangent;
else
    v=[];
end
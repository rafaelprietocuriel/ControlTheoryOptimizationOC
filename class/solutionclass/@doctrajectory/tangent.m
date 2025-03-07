function v=tangent(docTrj,varargin)


if isfield(docTrj.solverinfo,'tangent')
    v=docTrj.solverinfo.tangent;
else
    v=[];
end
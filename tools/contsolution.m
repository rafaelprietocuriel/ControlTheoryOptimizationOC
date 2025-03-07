function [contSol,contClass]=contsolution(contRes,varargin)
%
% CONT2IDX returns indices where a curve crosses a specific parameter value
%
% IDX=CONT2IDX(CONTPAR,PARVAL) returns the index/indices where the
% curve stored as a vector of discrete numbers CONTPAR, usually coming from
% a continuation process, crosses the parameter value PARVAL.

contSol=[];
contClass=[];
if ~isstruct(contRes) || ~isfield(contRes,'ContinuationSolution') || ~isfield(contRes,'ContinuationClassification')
    return
end

contSol=contRes.ContinuationSolution;
contClass=contRes.ContinuationClassification;
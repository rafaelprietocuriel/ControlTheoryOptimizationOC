function n=multiplicity(ocMultiPath)
%
%
if isempty(ocMultiPath)
    n=0;
else
    n=numel(ocMultiPath.solutionclass);
end
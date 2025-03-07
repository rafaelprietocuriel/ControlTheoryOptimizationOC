function n=multiplicity(docMultiPath)
%
%
if isempty(docMultiPath)
    n=0;
else
    n=numel(docMultiPath.solutionclass);
end
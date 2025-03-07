function out=timehorizon(ocMultiPath)

n=multiplicity(ocMultiPath);
out=zeros(1,n);
for ii=1:n
    out(ii)=timehorizon(ocMultiPath.solutionclass{ii});
end

function ocMultiPath=change(ocMultiPath,varargin)

n=multiplicity(ocMultiPath);
out=zeros(1,n);
for ii=1:n
    ocMultiPath.solutionclass{ii}=change(ocMultiPath.solutionclass{ii},varargin{:});
end
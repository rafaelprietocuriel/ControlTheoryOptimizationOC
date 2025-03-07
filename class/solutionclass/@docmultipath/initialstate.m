function out=initialstate(docMultiPath)

for ii=1:length(docMultiPath.solutionclass)
    out(:,ii)=initialtime(docMultiPath.solutionclass{ii});
end
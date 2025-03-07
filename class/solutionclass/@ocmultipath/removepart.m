function ocMP=removepart(ocMP,part)

counter=0;
for ii=1:multiplicity(ocMP)
    if ~ii==part
        counter=counter+1;
        ocMP.solutionclass{counter}=ocMP.solutionclass{ii};
    end
end

function ocMP=redefinearc(ocMP,newposition,arcid)

for ii=1:multiplicity(ocMP)
    ocMP.solutionclass{ii}=redefinearc(ocMP.solutionclass{ii},newposition{ii},arcid{ii});
end

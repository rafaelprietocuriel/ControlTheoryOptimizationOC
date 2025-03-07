function ocgTrj=real(ocgTrj)

for ii=1:length(ocgTrj.octrajectory.y)
    ocgTrj.octrajectory.y{ii}=real(ocgTrj.octrajectory.y{ii});
end
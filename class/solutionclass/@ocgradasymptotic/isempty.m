function b=isempty(ocgAsym)
%
%
if multiplicity(ocgAsym)>1
    b=1;
    for ii=1:multiplicity(ocgAsym)
        b=b&&isempty(ocgAsym(ii).ocgradtrajectory);
    end
else
    b=isempty(ocgAsym.ocgradtrajectory);
end

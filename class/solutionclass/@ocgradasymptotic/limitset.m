function out=limitset(ocgAsym)

out=ocgAsym.limitset;
for ii=2:multiplicity(ocgAsym)
    out(ii)=ocgAsym(ii).limitset;
end

function b=isequilibrium(ocgLim)
%
%
b=zeros(1,multiplicity(ocgLim));
if multiplicity(ocgLim)>1
    for ii=1:multiplicity(ocgLim)
        b(ii)=~ocgLim(ii).period;
    end
else
    b=~ocgLim.period;
end
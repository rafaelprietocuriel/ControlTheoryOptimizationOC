function b=isempty(ocgTrj)
%
%
if multiplicity(ocgTrj)>1
    b=1;
    for ii=1:multiplicity(ocgTrj)
        b=b&&isempty(ocgTrj(ii).modelname);
    end
else
    b=isempty(ocgTrj.modelname);
end

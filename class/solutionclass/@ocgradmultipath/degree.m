function out=degree(ocgMTrj)
out=[];
if isempty(ocgMTrj)
    return
end
out=ocgMTrj.multiplicity.number;
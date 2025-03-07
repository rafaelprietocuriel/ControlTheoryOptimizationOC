function out=multiplicity(ocgMTrj)

switch octype(ocgMTrj)
    case 'static'
        out=ocgMTrj.multiplicity.number;
    otherwise
        out=ocgMTrj.multiplicity;
end
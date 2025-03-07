function b=isempty(ocgMTrj)
%
%
switch octype(ocgMTrj)
    case 'static'
        b=isempty(multiplicity(ocgMTrj));
    otherwise
        b=~isfield(multiplicity(ocgMTrj),'number');
end

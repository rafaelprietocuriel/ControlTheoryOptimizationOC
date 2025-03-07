function term=ocmatsimple(term,symkernel)
% OCMATSIMPLE  

if nargin==1
    symkernel='maple';
end

switch symkernel
    case 'maple'
        term=maple('simplify',term);
    case 'mupad'
        term=char(simple(sym(term)));
end
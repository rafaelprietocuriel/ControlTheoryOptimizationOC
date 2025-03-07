function term=ocmatevalf(term,symkernel)

if nargin==1
    symkernel='maple';
end

switch symkernel
    case 'maple'
        term=maple('evalf',term);
end

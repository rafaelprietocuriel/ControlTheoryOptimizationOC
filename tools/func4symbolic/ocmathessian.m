function H=ocmathessian(f,x,symkernel)
%ocmatjacobian
H='';
if nargin==2
    symkernel='maple';
end

switch symkernel
    case 'maple'
        H=maple('hessian',f,x);
    case 'mupad'
        if verLessThan('symbolic','8')
            symx=sym(x);
            H=char(jacobian(jacobian(sym(f),symx),symx));
        else
            symx=mystr2sym(x);
            %H=char(jacobian(jacobian(str2sym(f),symx),symx));
            H=char(hessian(mystr2sym(f),symx));
        end
end


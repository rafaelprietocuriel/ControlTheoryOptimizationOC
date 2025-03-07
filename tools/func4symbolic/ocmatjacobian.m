function J=ocmatjacobian(f,x,symkernel,opt)
% OCMATJACOBIAN calculates the Jacobian using commands corresponding to the
% used kernel of the symbolic toolbox

J='';
if nargin==2
    symkernel='maple';
    opt=defaultocoptions;
end
if nargin==3
    opt=defaultocoptions;
end
Simplify=strcmp(getocoptions(opt,'INIT','Simplify'),'on');

switch symkernel
    case 'maple'
        J=maple('jacobian',f,x);
        if Simplify
            J=maple('simplify',J,'power');
        end
    case 'mupad'
        if verLessThan('symbolic','8')
            J=jacobian(sym(f),sym(x));
            if Simplify
                J=simple(J);
            end
            J=char(J);
        else
            J=jacobian(mystr2sym(f),mystr2sym(x));
            if Simplify
                J=simple(J);
            end
            J=strrep(char(J),';','],[');
        end
end


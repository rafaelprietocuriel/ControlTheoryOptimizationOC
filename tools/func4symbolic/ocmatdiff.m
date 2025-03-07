function J=ocmatdiff(f,x,symkernel)
% OCMATJACOBIAN calculates the Jacobian using commands corresponding to the
% used kernel of the symbolic toolbox

J='';
if nargin==2
    symkernel='maple';
end

switch symkernel
    case 'maple'
        if strncmp(f,'matrix',6)
            J=strrep(char(diff(f,x)),'],[','], [');
        else
            J=maple('diff',f,x);
        end
    case 'mupad'
        if verLessThan('symbolic','8')
            J=char(diff(sym(f),sym(x)));
        else
            J=char(diff(str2sym(f),str2sym(x)));
        end
end


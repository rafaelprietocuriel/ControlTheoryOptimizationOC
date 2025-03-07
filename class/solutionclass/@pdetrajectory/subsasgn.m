function pdeTrj=subsasgn(pdeTrj,index,rhs)
%
%
flag=0;
if length(index)==1
    switch index.type
        case '.'
            if strcmp(index.subs,'linearization')
                if size(rhs,3)==length(pdeTrj.x) && size(rhs,1)==size(rhs,2) && size(rhs,2)==size(pdeTrj.y,1)
                    pdeTrj.linearization=rhs;
                else
                    ocmatmsg('Linearization has the wrong size.')
                end
            else
                flag=1;
            end
        otherwise
            flag=1;
    end
        
else
    flag=1;
end

if flag
    ocmatmsg('Assigning is not implemented.')
end
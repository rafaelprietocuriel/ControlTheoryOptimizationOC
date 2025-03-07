function out=subsref(pdeTrj,index)
%
%
flag=0;
if length(index)==1
    switch index.type
        case '.'
            out=pdeTrj.(index.subs);
        otherwise
            flag=1;
    end
        
else
    flag=1;
end

if flag
    ocmatmsg('Referencing is not implemented.')
end

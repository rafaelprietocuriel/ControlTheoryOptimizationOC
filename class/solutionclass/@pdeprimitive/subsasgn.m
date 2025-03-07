function pdePrim=subsasgn(pdePrim,index,rhs)
%
%
flag=0;
if length(index)==1
    switch index.type
        case '.'
            if strcmp(index.subs,'linearization')
                pdePrim.pdetrajectory.linearization=rhs;
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
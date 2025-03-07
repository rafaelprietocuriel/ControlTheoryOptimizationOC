function out=subsref(pdeAsym,index)
%
%
flag=0;

out=pdeAsym;
for ii=1:length(index)
    actindex=index(ii);
    switch actindex.type
        case '.'
            if ii==1
                if strcmp(actindex.subs,'pdeprimitive')
                    out=out.pdeprimitive;
                elseif strcmp(actindex.subs,'pdetrajectory')
                    out=out.pdetrajectory;
                else
                    flag=1;
                end
            elseif ii==2
                out=out.(actindex.subs);
            end
        otherwise
            flag=1;
    end
end

if flag
    ocmatmsg('Referencing is not implemented.')
end

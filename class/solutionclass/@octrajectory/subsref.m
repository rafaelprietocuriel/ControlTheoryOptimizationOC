function out=subsref(ocTrj,index)
%
%
try
    if length(index)==1
        if strcmp(index.type,'()')
            if length(index.subs)==1
                out=ocTrj;
            else
                out=[];
            end
            return
        else
            out=subsref(struct(ocTrj),index);
        end
    else
        out=subsref(struct(ocTrj),index);
    end
catch
    [msgstr, msgid]=lasterr;
    rethrow(lasterror)
end
function out=subsref(ocComp,index)
%
%
switch numel(index)
    case 1
        if strcmp(index.type,'()')
            try
                if length(index.subs)==1
                    out=ocComp.path{index.subs{1}};
                else
                    out=ocComp.path{index.subs{1}}(index.subs{2});
                end
                return
            catch
                [msgstr, msgid]=lasterr;
                rethrow(lasterror)
            end
        end
end
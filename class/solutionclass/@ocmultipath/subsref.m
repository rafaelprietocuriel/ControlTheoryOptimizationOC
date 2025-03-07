function out=subsref(ocMP,index)
%
%
out=[];
switch numel(index)
    case 1
        if strcmp(index.type,'()')
            try
                out=ocMP.solutionclass{index.subs{1}};
                return
            catch
                [msgstr, msgid]=lasterr;
                rethrow(lasterror)
            end
        elseif strcmp(index.type,'.')
            try
                out=subsref(struct(ocMP),index);
                return
            catch
                [msgstr, msgid]=lasterr;
                rethrow(lasterror)
            end
        else
            ocmaterror('')
        end
    otherwise
        if strcmp(index(1).type,'()')
            try
                out=subsref(ocMP.solutionclass{index(1).subs{1}},index(2:end));
                return
            catch
                [msgstr, msgid]=lasterr;
                rethrow(lasterror)
            end
        else strcmp(index(1).type,'.')
            if strcmp(index(2).type,'{}')
                out=ocMP.(index(1).subs){index(2).subs{1}};
            end
        end
end
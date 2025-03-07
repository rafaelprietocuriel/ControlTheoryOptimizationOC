function out=subsref(docMultiPath,index)
%
%
out=[];
switch numel(index)
    case 1
        if strcmp(index.type,'()')
            try
                out=docMultiPath.solutionclass{index.subs{1}};
                return
            catch
                [msgstr, msgid]=lasterr;
                rethrow(lasterror)
            end
        elseif strcmp(index.type,'.')
            try
                out=subsref(struct(docMultiPath),index);
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
                out=subsref(docMultiPath.solutionclass{index(1).subs{1}},index(2:end));
                return
            catch
                [msgstr, msgid]=lasterr;
                rethrow(lasterror)
            end
        else
        end
end
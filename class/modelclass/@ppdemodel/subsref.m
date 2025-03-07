function out=subsref(ppdeObj,index)
%
%
try
    out=subsref(struct(ppdeObj),index);
catch
    [msgstr, msgid]=lasterr;
    if strcmp(msgid,'MATLAB:nonExistentField')
        ocmatmsg('Class ''%s'' only consists of fields ''Model'' and ''Result''!\n',class(ppdeObj))
    end
    if strcmp(msgid,'MATLAB:badsubscript')
        ocmatmsg('There is no array form of class ''%s''!\n',class(ppdeObj))
    end
    rethrow(lasterror)
end
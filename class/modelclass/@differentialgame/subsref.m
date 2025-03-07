function out=subsref(dgObj,index)
%
%
try
    out=subsref(struct(dgObj),index);
catch
    [msgstr, msgid]=lasterr;
    if strcmp(msgid,'MATLAB:nonExistentField')
        ocmatmsg('Class ''%s'' only consists of fields ''Model'' and ''Result''!\n',class(dgObj))
    end
    if strcmp(msgid,'MATLAB:badsubscript')
        ocmatmsg('There is no array form of class ''%s''!\n',class(dgObj))
    end
    rethrow(lasterror)
end
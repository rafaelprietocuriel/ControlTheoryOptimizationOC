function out=subsref(odeObj,index)
%
%
try
    out=subsref(struct(odeObj),index);
catch
    [msgstr, msgid]=lasterr;
    if strcmp(msgid,'MATLAB:nonExistentField')
        ocmatmsg('Class ''%s'' only consists of fields ''Model'' and ''Result''!\n',class(odeObj))
    end
    if strcmp(msgid,'MATLAB:badsubscript')
        ocmatmsg('There is no array form of class ''%s''!\n',class(odeObj))
    end
    rethrow(lasterror)
end
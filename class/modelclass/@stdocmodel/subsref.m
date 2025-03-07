function out=subsref(ocObj,index)
%
%
try
    out=subsref(struct(ocObj),index);
catch
    [msgstr, msgid]=lasterr;
    if strcmp(msgid,'MATLAB:nonExistentField')
        ocmatmsg('Class ''%s'' only consists of fields ''Model'' and ''Result''!\n',class(ocObj))
    end
    if strcmp(msgid,'MATLAB:badsubscript')
        ocmatmsg('There is no array form of class ''%s''!\n',class(ocObj))
    end
    rethrow(lasterror)
end
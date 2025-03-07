function out=subsref(ocCuv,index)
%
%
try
    out=subsref(struct(ocCuv),index);
catch
    [msgstr, msgid]=lasterr;
    rethrow(lasterror)
end
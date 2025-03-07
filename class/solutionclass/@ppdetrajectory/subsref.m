function out=subsref(ppdeTrj,index)
%
%
try
    out=subsref(struct(ppdeTrj),index);
catch
    [msgstr, msgid]=lasterr;
    rethrow(lasterror)
end
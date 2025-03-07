function out=subsref(ocTrj,index)
%
%
try
    out=subsref(struct(ocTrj),index);
catch
    [msgstr, msgid]=lasterr;
    rethrow(lasterror)
end
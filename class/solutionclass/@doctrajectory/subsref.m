function out=subsref(docTrj,index)
%
%
try
    out=subsref(struct(docTrj),index);
catch
    [msgstr, msgid]=lasterr;
    rethrow(lasterror)
end
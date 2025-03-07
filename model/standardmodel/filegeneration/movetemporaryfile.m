function [status,message,messageid]=movetemporaryfile(temporaryfilename,newfilwname)

[status,message,messageid]=movefile(temporaryfilename,newfilwname,'f');
if ~status
    try
        fclose('all')
    end
    [status,message,messageid]=movefile(temporaryfilename,newfilwname,'f');
end
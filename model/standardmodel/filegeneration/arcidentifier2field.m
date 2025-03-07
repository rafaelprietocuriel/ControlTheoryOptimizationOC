function arcfield=arcidentifier2field(arcidentifier)
arcfield=[];
if ischar(arcidentifier)
    if strncmp(arcidentifier,'arc',3)
        arcfield=arcidentifier;
        return
    end
    arcfield=['arc' arcidentifier];
    return
end

if iscell(arcidentifier)
    arcfield=['arc' arcidentifier{1}];
    for ii=2:numel(arcidentifier)
        arcfield=[arcfield 'p' arcidentifier{ii}];
    end
end
function arcidentifier=field2arcidentifier(arcfield)
arcidentifier=[];
if isempty(arcfield)
    return
end
arcfield=strrep(arcfield,'arc','');

arcidentifier=regexp(arcfield,'p','split');

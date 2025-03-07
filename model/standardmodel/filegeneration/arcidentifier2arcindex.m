function arcindex=arcidentifier2arcindex(arcidentifier)
arcindex=[];
if iscell(arcidentifier)
    arcindex=str2num(char(arcidentifier));
    arcindex=arcindex(:).';
elseif ischar(arcidentifier)
    arcindex=str2double(arcidentifier);
end
arcindex=arcindex+1;
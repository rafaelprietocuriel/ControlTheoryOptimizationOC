function arcindex=arcarg2arcindex(arcarg)
%
if ~isnumeric(arcarg)
    ocmaterror('Input argument has to be a non negative integer array.')
end
arcindex=arcarg+1;
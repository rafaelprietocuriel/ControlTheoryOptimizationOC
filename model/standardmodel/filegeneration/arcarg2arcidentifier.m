function arcid=arcarg2arcidentifier(arcarg)
%
if ~isnumeric(arcarg)
    ocmaterror('Input argument has to be a non negative integer array.')
end
arcid=num2str(arcarg);
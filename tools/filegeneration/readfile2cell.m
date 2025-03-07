function filetext=readfile2cell(filename2read,scanformat,varargin)
%
% READFILE2CELL reads the content of an (ASCII) file and returns each line
% in a separate cell. This command is used for processing the
% initialization file and for the generation of model files from template 
% files.  

fid_read=fopen(filename2read,'r');
if fid_read==-1
    ocmaterror(['Problems to open file: ''%s'''],filename2read)
end
filetext=textscan(fid_read,scanformat,varargin{:});
filetext=filetext{1};
fclose(fid_read);

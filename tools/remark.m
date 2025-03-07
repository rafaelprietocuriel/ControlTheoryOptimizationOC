function varargout=remark(basicfilename,folder,varargin)
%
% REMARK creates or opens a text file for remarks.
%
% REMARK(FNAME) the name of text file is given by the basic filename FNAME
% appending Remark and the extension 'ocr'. This file is created or opened
% if it already exists in the current directory. The actual date is
% appended at the end of the file. On a 'PC' it is opened with a user
% specified editor, if the extension is known to the system, e.g,
% Notepad++, or the MATLAB editor in any other cases. This command shall
% help the user to make model specific notes during her/his analysis.
%
% REMARK(FNAME,FOLDER) with FOLDER a specific folder for the file can be
% specified.
%
% REMARK(FNAME,FOLDER,TYPE) the TYPE specifies the extension of the file.
%   TYPE 
%       'standardmodel'         : 'ocr' 
%       'uncontrolledodemodel'  : 'odr' 
%       'spacedistributedmodel' : 'pdr' 
% and 'txt' otherwise.        
%
% REMARK(FNAME,FOLDER,TYPE,WRITEDATE) the WRITEDATE specifies if a string
% with the actual is appended. WRITEDATE = 0 or 1.
%
% FIDWRITE=REMARK(...) the file identifier FIDWRITE is returned.



basicextension='';
writedate=[];
if nargin==1
    folder='./';
end
if nargin>=3
    basicextension=varargin{1};
end
if nargin>=4
    writedate=varargin{2};
end
if isempty(basicextension)
    extension='txt';
else
    extension=[basicextension 'r'];
end
if isempty(writedate)
    writedate=true;
end
textfile=fullfile(folder,[basicfilename 'Remark.' extension]);
if ~exist(textfile,'file')
    fprintf('File ''%s'' does not exist.\n',textfile)
    answer=input('Create it? (y)/n: ','s');
    while 1
        if isempty(answer)
            % default value 'y'
            answer='y';
        end
        if strcmpi(answer,'n')
            if nargout==1
                varargout{1}=[];
            end
            return
        elseif strcmpi(answer,'y')
            break
        end
    end

end
if writedate || ~exist(textfile,'file')
    [fidwrite,msg]=fopen(textfile,'a+');
    if fidwrite==-1
        ocmatmsg(msg)
        return
    end
    fprintf(fidwrite, '\n%s\n',repmat('%',1,75));
    fprintf(fidwrite, '%s\t%s:\n\n',repmat('%',1,2),datestr(date,'dd-mmm-yyyy'));
    fclose(fidwrite);
else
    fidwrite=[];
end
if ispc
    try
        winopen(textfile);
    catch
        edit(textfile)
    end
else
    edit(textfile)
end
if nargout==1
    varargout{1}=fidwrite;
end
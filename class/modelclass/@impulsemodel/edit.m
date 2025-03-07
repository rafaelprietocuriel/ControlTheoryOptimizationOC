function edit(ocObj,varargin)
%
% EDIT

if nargin==1
    what='init';
else
    what=varargin{1};
end

switch lower(what(1))
    case 'i' % initialization file
        textfile=fullocmatfile(getocmatfolder('initializationfile'),initfilename(ocObj));
        if ~exist(textfile,'file')
            ocmatmsg('Initialization file is not in the default directory. Search in MATLAB path.\n')
            textfile=initfilename(ocObj);
            if ~exist(textfile,'file')
                ocmatmsg('Initialization file does not exist on MATLAB path.\n')
                return
            end
        end
    case 'r' % remark file
        remarkfile=[modelname(ocObj) 'Remark.' basicextension(ocObj) 'r'];
        textfile=fullocmatfile(getocmatfolder(ocObj,'userdata'),remarkfile);
        if ~exist(textfile,'file')
            answer=input([remarkfile ' is not on MATLAB path or does not exist. Create it?  (y)/n: '],'s');
            if isempty(answer)
                % default value 'y'
                answer='y';
            end
            if strcmpi(answer,'y')
                remark(ocObj)
            end
            textfile='';
        end
        
    case 'h'
        textfile=fullocmatfile(userdatafolder(ocObj),[modelname(ocObj) 'History.m']);
        if ~exist(textfile,'file')
            answer=input(['Model history file is not on MATLAB path or does not exist. Create it?  (y)/n: '],'s');
            if isempty(answer)
                % default value 'y'
                answer='y';
            end
            if strcmpi(answer,'y')
                history(ocObj,'open')
            end
            textfile='';
        end
    otherwise
        textfile='';
end
if ~isempty(textfile)
    if ispc
        try
            winopen(textfile);
        catch
            edit(textfile)
        end
    else
        edit(textfile)
    end
end

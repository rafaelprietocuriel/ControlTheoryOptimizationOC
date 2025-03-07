function folder=getocmatfolder(foldertype,modeltype,modelname)
% returns a folder within the ocmat folder structure
%
% folder=GETOCMATFOLDER()
folder='';
if nargin<2
    modeltype='standardmodel';
end
if nargin<3
    modelname='';
end


switch foldertype
    case 'model'
        folder='model';

    case 'out'
        folder=fullfile(getocmatfolder('usermodel',[]),'out');
        
    case 'initializationfile'
        folder=fullfile(getocmatfolder('model',[]),'initfiles');

    case 'templatefolder'
        folder=fullfile(getocmatfolder('model',[]),modeltype,'templatefiles');
        
    case 'usermodel'
        folder=fullfile(getocmatfolder('model',[]),'usermodel');

    case 'specificmodel'
        if isempty(modelname)
            ocmatmsg('No modelname provided.')
            return
        end
        folder=fullfile(getocmatfolder('usermodel',[]),modelname);
        
    case 'userdata'
        folder=fullfile(getocmatfolder('specificmodel',[],modelname),'data');

    case 'latexdoc'
        folder=fullfile('doc','latex');

    case 'matlabdoc'
        folder=fullfile(getocmatfolder('latexdoc',[]),'matlab');
end
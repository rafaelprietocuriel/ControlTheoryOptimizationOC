function [zipfile,status]=compressmodel(ocObj,varargin)
%

savedir=fullfile(fullocmatfile(modelfolder(ocObj)),'*');
zipfile=[fullfile(fullocmatfile(getocmatfolder('out')),modelname(ocObj)) '-' datestr(now, 'dd-mm-yyyy') '.zip'];
initializationfile=fullfile(fullocmatfile(getocmatfolder('initializationfile')),initfilename(ocObj));
status=system(['zip -r ' zipfile ' ' savedir ' ' initializationfile]);

function [status,result]=copyocmat(varargin)

if strcmp(computer,'PCWIN') || strcmp(computer,'PCWIN64')
    actdir=cd;
    ocmatpath=getocmatpath();
    cd(ocmatpath)
    copyocmatbatchfile=fullfile(ocmatpath,'copyocmat.bat');
    for ii=1:nargin
        copyocmatbatchfile=[copyocmatbatchfile ' ' varargin{ii}];
    end
    [status,result]=dos(copyocmatbatchfile,'-echo');
    cd(actdir)
else
end
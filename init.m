function init(varargin)
% INIT set the path variable for ocmat
%
% INIT(modelname1,...) includes the model(s) modelname1,...
%

actualpath=cd;
folder={'class','tools','options','cont',fullfile('model','standardmodel'),fullfile('model','default')};
for ii=1:length(folder)
    addpath(genpath(fullfile(actualpath,folder{ii})))
end
addpath(fullfile(actualpath,getocmatfolder('usermodel'),'out'))
for ii=1:nargin
    addpath(genpath(fullfile(actualpath,getocmatfolder('specificmodel',[],varargin{ii}))))
end
savepath
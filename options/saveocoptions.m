function saveocoptions(opts,filename,varargin)
%
% SAVEOCOPTIONS saves a certain optionstructure.
%
% SAVEOCOPTIONS(OPTS,FILENAME) The option structure OPTS will be saved in the
% file called FILENAME. The default filepath is data/options.
%
% SAVEOCOPTIONS(OPTS,FILENAME,PATH) The option structure OPTS will be saved in
% file called FILENAME, located at the PATH defined in the input arguments. 

if ~isstruct(opts)
    error('MATLAB: saveocoptions: The options must be within a structure.');
end

if ~isstr(filename)
    error('MATLAB: saveocoptions: The filename must be a string.');
end

% Evaluate whether default or individual path is supposed to be used.

if nargin==3
    path=varargin{1};
    if isstr(path)
        file=[path filesep filename];
    else
        error('MATLAB: saveocoptions: The filepath must be a string.');
    end
else
    path=getocpath('data');
    file=[path 'options' filesep filename];
end

save(file,'opts');

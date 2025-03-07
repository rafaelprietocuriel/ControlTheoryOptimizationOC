function opts=loadocoptions(filename,varargin)
%
% LOADOCOPTIONS An option structure that has priorily been saved to the
% file filename will be loaded.
%
% OPTS=LOADOCOPTIONS(FILENAME) The option structure contained in FILENAME,
% located at the defaultpath at data/options/ will be returned in the
% option structure OPTS.  
%
% OPTS=LOADOCOPTIONS(FILENAME,PATH) The option structure contained in
% FILENAME, located at the PATH defined in the input arguments will be
% returned in the option structure OPTS. 


if ~ischar(filename)
    error('MATLAB: loadocoptions: filename must be a string.')
end

% Evaluate whether the required file where this option is saved is to be
% found at the defaultpath or somewhere else.

if nargin==2
    file=[varargin{1} filesep filename];
else
    path=getocpath('data');
    file=[path 'options' filesep filename];
end

op=load(file);  
opts=op.opts;

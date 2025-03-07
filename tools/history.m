function history(modelname,varargin)

savedir='';
fn='';

if nargin>=2
   savedir=varargin{1};
end

if nargin>=3
   fn=varargin{2};
end

if isempty(fn)
    fn=[modelname 'History.m'];
end
if isempty(savedir)
    savedir=fullfile(getocmatpath,'history'); % Data directory
end
if verLessThan('Matlab','8.3.0')
    hfile=fullfile(prefdir,'history.m');
else
    hfile=fullfile(prefdir,'History.xml');
end
if ~exist(hfile,'file')
    ocmatmsg('Cannot open history file.\n')
    return
end
mathist = fileread(hfile);
if ~verLessThan('Matlab','R2014b')
    mathist = regexprep(mathist, '(<[^>]+>\s*)+', '\n', 'lineanchors');
    % translate html entities and remove leading newline
    mathist = strrep(mathist(3:end), '&gt;', '>');
    mathist = strrep(mathist, '&lt;', '<');
end
fidwrite = fopen(fullfile(savedir,fn),'a+');
if fidwrite==-1
    ocmatmsg('Cannot open file to write %s.\n',fullfile(savedir,fn))
    return
end

fwrite(fidwrite,mathist);
fclose(fidwrite);



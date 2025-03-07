function showocoptions(varargin)
%
% SHOWOCOPTIONS is a function to display the available categories, names and
% values of an option structure.
%
% SHOWOCOPTIONS shows the defaultoptions.
%
% SHOWOCOPTIONS(OPTS) displays the options of the option structure OPTS.
% 
% SHOWOCOPTIONS(CAT) displays the content of the option category CAT of
% the default options
%
% SHOWOCOPTIONS(OPTS,CAT) displays the content of the option category CAT
% of the option structure OPTS.


switch nargin
    case 0
       opts=defaultocoptions();
       showalloptions(opts);
    case 1
        if isstruct(varargin{1})
            opts=varargin{1};
            showalloptions(opts);
        elseif isfield(defaultocoptions,varargin{1})
            dispcat(varargin{1},getfield(defaultocoptions,varargin{1}));
        else
            error('Invalid input argument');
        end
    case 2
        if isstruct(varargin{1}) && isfield(varargin{1},varargin{2})
            dispcat(varargin{2},getfield(varargin{1},varargin{2}));
        else
            error('Invalid input argument');
        end
    otherwise
        error('Too many input arguments.')
end
end
            
            
function showalloptions(opts)

disp('---------------------------OC Options-----------------------');

fname=fields(opts);

for ii=1:length(fname)
    dispcat(char(fname(ii)),getfield(opts,char(fname(ii))));
end
disp('--------------------------------------------------');
end

function dispcat(catname,catcontent)

underl='';

    for uu=1:length(catname)
        underl=strcat(underl,'~');
    end

    % Displays category names
    disp(upper(catname));
    disp(underl);
    disp(catcontent);
    disp(' ');
end

        

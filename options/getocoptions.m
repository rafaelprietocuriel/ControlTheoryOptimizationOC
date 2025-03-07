function varargout = getocoptions(varargin)
%
% GETOCOPTIONS returns the value of the options defined in varargin
%
% [VALUE1 VALUE2 ...]=GETOCOPTIONS(CAT1,NAME1,CAT2,NAME2,...) returns the
% default values of the options with the named NAME1, NAME2,... located in
% the categories CAT1,CAT2,.. 
%
% [VALUE1 VALUE2 ...]=GETOCOPTIONS(CAT1,NAME1,DEFAULTVALUE1,CAT2,NAME2,...)
% returns the default values of the options with the named NAME1, NAME2,...
% located in the categories CAT1,CAT2,.. . If the default values are empty
% the optional input argument defaultvalue assigns an alternative value to
% this option. This default value does not have to be assigned for every
% option that is to be returned.
%
% [VALUE1 VALUE2 ...]=GETOCOPTIONS(OPT,CAT1,NAME1,CAT2,NAME2,...) returns the
% value of the options with the names NAME1,NAME2,... of the categories
% CAT1,CAT2,... of the option structure OPT.
%
% [VALUE1 VALUE2 ...]=GETOCOPTIONS(OPTS,CAT1,NAME1,DEFAULTVALUE1,CAT2,NAME2,...)
% returns the value of the options with the names NAME1,NAME2,... of the
% categories CAT1,CAT2,... of the option structure OPT. To avoid an empty
% output argument one can assign default values that are returned in case
% the required value of the option is empty in this option structure.
%
% [VALUE1 VALUE2 ...]=GETOCOPTIONS(OPTS,CAT1,NAME1,DEFAULTVALUE1,NAME2,CAT2,...)
% returns the value of the options with the names name1,name2,... of the
% categorie CAT1 of the option structure OPT. To avoid an empty
% output argument one can assign default values that are returned in case 
% the required value of the option is empty in this option structure.

optflag=2;

% Evaluation which option structure is to be accessed.

if isstruct(varargin{1})
    opts=varargin{1};
elseif isempty(varargin{1})
    opts=defaultocoptions;
else
    opts=defaultocoptions;
    optflag=1;
end


% Evaluation of which options are required.

idx=1;
cat=[];
catcont=[];
optcont=[];
ii=optflag;

while ii<=nargin
    % Check whether input argument is an option category
    if isfield(opts,varargin{ii})
        cat=varargin{ii};
        catcont=getfield(opts,cat);
        ii=ii+1;
    else
        if isempty(cat)
            error('Category name not set.')
        end
        
        if isfield(catcont,varargin{ii})
            optcont=getfield(catcont,varargin{ii});
            
            if ii+1 <= nargin &&  ~isfield(opts,varargin{ii+1}) && ~isfield(catcont,varargin{ii+1})&& ~isfield(optcont,varargin{ii+1})                  
               if isempty(optcont)
                  optcont=varargin{ii+1};
               end                   
               ii=ii+2;
            elseif ~isstruct(optcont)
                ii=ii+1;
            end
            
            bb=0;
            subopt=[];
            if isstruct(optcont)              
                nel=1;
                while nel
                    
                    if (ii+1<=nargin && ~isfield(optcont,varargin{ii+1})&&isfield(opts,varargin{ii+1}))
                        nel=0;
                        ii=ii+1;
                    elseif ii+1 <=nargin && isfield(optcont,varargin{ii+1})
                        subopt=getfield(optcont,varargin{ii+1});
                        if isempty(subopt)
                            bb=1;
                        end
                        if ii+2 <= nargin &&  ~isfield(opts,varargin{ii+2}) && ~isfield(catcont,varargin{ii+2}) && ~isfield(optcont,varargin{ii+2})                  
                            if isempty(subopt)
                                subopt=varargin{ii+2};
                            end                   
                            ii=ii+2;
                        else
                            ii=ii+1;
                        end
                        varargout{idx}=subopt;
                        idx=idx+1; 
                    else
                        ii=ii+1;
                        nel=0;
                    end
                end
            end
            if isempty(subopt)&&isstruct(optcont)
                ii=ii+1;
            end
             if isempty(subopt) && ~bb                    
                varargout{idx}=optcont;
                idx=idx+1;
             end
        else
            error('Invalid option name');
        end
                       
    end        
end
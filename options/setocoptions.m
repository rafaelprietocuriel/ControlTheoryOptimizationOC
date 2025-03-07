function opts=setocoptions(varargin)
%
% SETOCOPTIONS lets you adjust the options of the OCMAT toolbox.
%
% OPT=SETOCOPTIONS(CAT1,NAME1,VALUE1,CAT2,NAME2,VALUE2,...) replaces the
% values of the options NAME1, NAME2,... in the categories CAT1,CAT2,.. of
% the default option structure with new values VALUE1, VALUE2,... and
% returns an option structure OPT with the changed options.
%
% OPT=SETOCOPTIONS(OPS,CAT1,NAME1,VALUE1,CAT2,NAME2,VALUE2,...) replaces
% the values of the options NAME1, NAME2,... in the categories CAT1,CAT2,..
% of the option structure OPT with new values VALUE1, VALUE2,... and
% returns the option structure OPT with the changed options. 
%
% OPT=SETOCOPTIONS(OPT,CAT1,NAME1,VALUE1,NAME2,VALUE2,CAT2,NAME3,VALUE3,..)
% replaces the values of the options NAME1, NAME2, ... of the category CAT1
% of the (optional) option structure OPT with the values VALUE1 and VALUE2
% and option NAME3 of category CAT2 with value VALUE3 etc. and returns
% the changed option structure OPT.

flag=1;

if isstruct(varargin{1})
    opts=varargin{1};
    flag=2;
elseif isempty(varargin{1})
    opts=defaultocoptions;
    flag=2;
else
    opts=defaultocoptions;
end

cat=[];
catcont=[];
ii=flag;

while ii<=nargin
    if isfield(opts,varargin{ii})
        cat=varargin{ii};
        catcont=getfield(opts,cat);
        ii=ii+1;
    else
        if isempty(cat)
            error('No category given.')
        end

        if isfield(catcont,varargin{ii})
            if ii+1<=nargin
                fname=varargin{ii};
                optcont=getfield(catcont,fname);
                nel=1;
                while nel
                    if isfield(optcont,varargin{ii+1})
                        if ii+2<=nargin
                            varargin{ii+1};
                            optcont=setfield(optcont,varargin{ii+1},varargin{ii+2});
                            catcont=setfield(catcont,fname,optcont);
                            opts=setfield(opts,cat,catcont);
                            if ii+3<=nargin && isfield(optcont,varargin{ii+3})
                                ii=ii+2;
                            else
                                ii=ii+3;
                                nel=0;
                            end
                        end
                    else
                        catcont=setfield(catcont,varargin{ii},varargin{ii+1});
                        opts=setfield(opts,cat,catcont);
                        ii=ii+2;
                        nel=0;
                    end
                end

            else
                error('Invalid number of input arguments.')
            end
        else
            error('Invalid option name.')
        end
    end

end
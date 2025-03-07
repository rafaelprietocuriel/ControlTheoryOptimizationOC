function out=mystr2sym(varargin)

if verLessThan('symbolic','8')
    if iscell(varargin{:})
        out=cell(nargin,1);
        for ii=1:nargin
            if nargin==1
                for jj=1:length(varargin{ii})
                    out{jj}=sym(varargin{ii}{jj});
                end
            end
        end
    else
        out=sym(varargin{:});
    end
else
    out=str2sym(['[' varargin{:} ']']);
end



function b=testsymbolic(S,symkernel)
%
% test if input argument S is a syntactically correct symbolic expression.
% S is a character of cell of characters

b=[];
if ischar(S)
    S=cellstr(S);
end
if ~iscell(S)
    ocmaterror('Input argument is not a character or cell.')
end
switch symkernel
    case 'maple'
        b=zeros(1,numel(S));
        for ii=1:numel(S)
            try
                if ~isempty(S{ii})
                    maple('evalf',S{ii});
                end
                b(ii)=1;
            end
        end
    case 'mupad'
        if verLessThan('symbolic','8')
            b=zeros(1,numel(S));
            for ii=1:numel(S)
                try
                    if ~isempty(S{ii})
                        mystr2sym(S{ii});
                    end
                    b(ii)=1;
                end
            end
        else
            b=zeros(1,numel(S));
            for ii=1:numel(S)
                try
                    if ~isempty(S{ii})
                        mystr2sym(S{ii});
                    end
                    b(ii)=1;
                end
            end
        end
    otherwise
        ocmatmsg('No kernel for the symbolic toolbox specified.\n')
end
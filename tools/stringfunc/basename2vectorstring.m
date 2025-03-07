function s=basename2vectorstring(c,idx,vectortype,varargin)
%
% MAKEVECTORSTRING composes string for symbolic calculations
%
% S=MAKESTR(C,N) composes string, build from character C and

s='';
separator=[];
if nargin>=4
    separator=varargin{1};
end
if isempty(idx)
    return
end
if isempty(separator)
    separator=',';
end
if isempty(separator)
    separator=',';
end
if ~ischar(c)
    ocmaterror('First argument is not a character')
end
lnsep=length(separator);
switch vectortype
    case 'index'
        for ii=idx
            s=[s c num2str(ii) separator];
        end
    case 'coord'
        for ii=idx
            s=[s c '(' num2str(ii) ')' separator];
        end
end
s(end-lnsep+1:end)=[];

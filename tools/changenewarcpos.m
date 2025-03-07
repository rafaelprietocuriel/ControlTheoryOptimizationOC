function newarcpos=changenewarcpos(newarcpos,d,n)

if nargin==1
    d='b';
    n=1;
end
if nargin==2
    n=1;
end

switch d
    case 'b'
        newarcpos=newarcpos+[-n;n];
    case 'l'
        newarcpos=newarcpos+[-n;0];
    case 'r'
        newarcpos=newarcpos+[0;n];
end
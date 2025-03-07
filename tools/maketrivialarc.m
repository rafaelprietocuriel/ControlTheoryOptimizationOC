function tarc=maketrivialarc(pt,arcarg,n)

if nargin==2
    n=5;
end
pt=pt(:);
tarc.x=linspace(0,1,n);
tarc.y=pt(:,ones(1,n));
tarc.arcinterval=[0 0];
tarc.arcarg=arcarg;
tarc.arcposition=[1 n]';
function out=arcargument(ppdeAsym,flag)

if nargin==1
    flag='ppdetrajectory';
end

switch flag
    case 'ppdetrajectory'
        out=arcargument(ppdeAsym.ppdetrajectory);
    case 'limitset'
        out=arcargument(ppdeAsym.limitset);
end

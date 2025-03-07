function out=timeinterval(ppdeAsym,flag)

if nargin==1
    flag='ppdetrajectory';
end

switch flag
    case 'ppdetrajectory'
        out=timeinterval(ppdeAsym.ppdetrajectory);
    case 'limitset'
        out=timeinterval(ppdeAsym.limitset);
end

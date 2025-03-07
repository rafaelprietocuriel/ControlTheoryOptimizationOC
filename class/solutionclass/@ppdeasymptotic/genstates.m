function out=genstates(ppdeAsym,flag)

if nargin==1
    flag='ppdetrajectory';
end

switch flag
    case 'ppdetrajectory'
        out=genstates(ppdeAsym.ppdetrajectory);
    case 'limitset'
        out=genstates(ppdeAsym.limitset);
end

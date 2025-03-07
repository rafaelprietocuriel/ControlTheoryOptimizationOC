function out=initialtime(ppdeAsym,flag)

if nargin==1
    flag='ppdetrajectory';
end

switch flag
    case 'ppdetrajectory'
        out=initialtime(ppdeAsym.ppdetrajectory);
    case 'limitset'
        out=initialtime(ppdeAsym.limitset);
end

function out=arcargument(ocAsym,flag)

if nargin==1
    flag='octrajectory';
end

switch flag
    case 'octrajectory'
        out=arcargument(ocAsym.octrajectory);
    case 'limitset'
        out=arcargument(ocAsym.limitset);
end

function out=arcposition(ocAsym,flag)

if nargin==1
    flag='octrajectory';
end

switch flag
    case 'octrajectory'
        out=arcposition(ocAsym.octrajectory);
    case 'limitset'
        out=arcposition(ocAsym.limitset);
end

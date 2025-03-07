function out=arcinterval(ocAsym,flag)

if nargin==1
    flag='octrajectory';
end

switch flag
    case 'octrajectory'
        out=arcinterval(ocAsym.octrajectory);
end

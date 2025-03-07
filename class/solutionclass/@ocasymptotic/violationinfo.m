function out=violationinfo(ocAsym,flag)

if nargin==1
    flag='octrajectory';
end

switch flag
    case 'octrajectory'
        out=violationinfo(ocAsym.octrajectory);
    case 'limitset'
        out=violationinfo(ocAsym.limitset);
end

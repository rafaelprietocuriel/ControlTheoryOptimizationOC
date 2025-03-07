function out=dependentvar(ocAsym,flag)

if nargin==1
    flag='octrajectory';
end

switch flag
    case 'octrajectory'
        out=dependentvar(ocAsym.octrajectory);
    case 'limitset'
        out=dependentvar(ocAsym.limitset);
end
function out=dependentvar(ocgAsym,flag)

if nargin==1
    flag='ocgtrajectory';
end

switch flag
    case 'ocgtrajectory'
        out=dependentvar(ocgAsym.ocgtrajectory);
    case 'limitset'
        out=dependentvar(ocgAsym.limitset);
end
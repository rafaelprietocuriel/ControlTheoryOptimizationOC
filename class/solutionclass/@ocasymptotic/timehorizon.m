function out=timehorizon(ocAsym,flag)

if nargin==1
    flag='octrajectory';
end

switch flag
    case 'octrajectory'
        out=timehorizon(ocAsym.octrajectory);
end

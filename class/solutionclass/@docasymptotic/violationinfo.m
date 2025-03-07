function out=violationinfo(docAsym,flag)

if nargin==1
    flag='doctrajectory';
end

switch flag
    case 'doctrajectory'
        out=violationinfo(docAsym.doctrajectory);
    case 'limitset'
        out=violationinfo(docAsym.limitset);
end

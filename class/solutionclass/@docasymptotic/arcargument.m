function out=arcargument(docAsym,flag)

if nargin==1
    flag='doctrajectory';
end

switch flag
    case 'doctrajectory'
        out=arcargument(docAsym.doctrajectory);
    case 'limitset'
        out=arcargument(docAsym.limitset);
end

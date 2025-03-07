function out=dependentvar(docAsym,flag)

if nargin==1
    flag='doctrajectory';
end

switch flag
    case 'doctrajectory'
        out=dependentvar(docAsym.doctrajectory);
        if strcmp(pathtype(docAsym),'u')
            out=out(:,end:-1:1);
        end
    case 'limitset'
        out=dependentvar(docAsym.limitset);
end
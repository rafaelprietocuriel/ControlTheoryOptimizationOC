function out=arcposition(docAsym,flag)

if nargin==1
    flag='doctrajectory';
end

switch flag
    case 'doctrajectory'
        out=arcposition(docAsym.doctrajectory);
        if strcmp(pathtype(docAsym),'u')
            out=out(2,end)-out([2 1],:)+1;
        end        
    case 'limitset'
        out=arcposition(docAsym.limitset);
end

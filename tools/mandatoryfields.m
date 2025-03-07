function out=mandatoryfields(classname)

switch classname
    case 'pdetrajectory'
        out={'x','y','femdata','arcarg','arcposition','arcinterval'};
    case 'pdeprimitive'
        out={'y','femdata','arcarg'};
    otherwise
        out='';
end
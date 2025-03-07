function out=dependentvar(ocTrj,varargin)

switch octype(solOb)
    case 'concentrated'
        out=ocTrj.variable.y;
end

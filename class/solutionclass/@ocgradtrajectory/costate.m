function out=costate(ocgTrj,varargin)
out=[];
if isempty(ocgTrj)
    return
end
if nargin==1
    octpe=octype(ocgTrj);
elseif nargin>1
    octpe=varargin{1};
end
switch octpe
    case 'concentrated'
        out=ocgTrj.variable.cst_y;
end
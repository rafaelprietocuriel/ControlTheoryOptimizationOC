function out=state(ocgMTrj,varargin)
out=[];
if isempty(ocgMTrj)
    return
end
if nargin==1
    octpe=octype(ocgMTrj);
elseif nargin>1
    octpe=varargin{1};
end
switch octpe
    case 'static'
        out=ocgMTrj.variable.y;
    case 'concentrated'
        out=ocgMTrj.variable.y;
        %degr=multiplicity(ocgMTrj);
        %out=reshape(out.',size(out,2)/degr.number,[],degr.number);
        %out=reshape(out,[],size(out,1),degr.number);
end
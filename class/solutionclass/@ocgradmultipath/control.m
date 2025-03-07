function out=control(ocgMTrj,varargin)
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
    case 'concentrated'
        out=ocgMTrj.variable.v;
        %degr=multiplicity(ocgMTrj);
        %out=reshape(out,[],size(out,2)/degr.number,degr.number);
end
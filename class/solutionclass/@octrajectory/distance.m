function ocTrj=distance(ocTrj1,ocTrj2,varargin)

spec=[];
coord=[];
if nargin>=3
    spec=varargin{1};
end
if nargin>=4
    coord=varargin{2};
end
if isempty(spec)
    spec=1;
end
if isempty(coord)
    coord=1:size(ocTrj1.y,1);
end
switch spec
    case 0
        ocTrj=ocTrj1;
        ocTrj.solverinfo=[];
        ocTrj.y=ocTrj1.y(coord,:)-ocTrj2.y(coord,:);
    case inf
    otherwise
        ocTrj=ocTrj1;
        ocTrj.solverinfo=[];
        ocTrj.y=sum(abs(ocTrj1.y(coord,:)-ocTrj2.y(coord,:)).^spec).^(1/spec);
end

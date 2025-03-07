function out=dependentvar(ocTrj,varargin)

pos=[];

if nargin>=2
    pos=varargin{1};
end
if isempty(pos)
    out=ocTrj.y;
else
    out=ocTrj.y(:,pos);
end
function [b,adminfo]=isadmissible(gdynPrim,ocObj,opt,userf,varargin)
%

if nargin==2
    opt=[];
    userf=[];
    equationfile=[];
end

if nargin==3
    userf=[];
end
if nargin>3
    equationfileidx=find(strcmpi(varargin,'equationfile'));
    if ~isempty(equationfileidx)
        equationfile=varargin{equationfileidx+1};
        varargin(equationfileidx+(0:1))=[];
    else
        equationfile='EquilibriumEquation';
    end
end
if isempty(equationfile)
    equationfile='EquilibriumEquation';
end
[b,adminfo]=isadmissible(dynprimitive(gdynPrim),ocObj,opt,userf,'equationfile',equationfile,varargin{:});
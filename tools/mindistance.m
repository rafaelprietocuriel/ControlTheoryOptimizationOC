function idx=mindistance(pts,refpt,varargin)
%

tol=[];
if nargin==3
    if isstruct(varargin{1})
        tol=getocoptions(varargin{1},'GENERAL','ZeroDeviationTolerance');
    else
        tol=varargin{1};
    end
end

if isempty(tol)
    tol=1e-5;
end
reltol=1e-5;
    
relval=(pts-refpt(:,ones(1,size(pts,2))))./(abs(pts)+reltol);
d=sqrt(sum(relval.^2));
idx=find(d<tol);

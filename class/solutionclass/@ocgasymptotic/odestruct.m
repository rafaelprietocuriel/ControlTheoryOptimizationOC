function sol=odestruct(ocgAsym,parflag)
%
%

sol=[];
if isempty(ocgAsym)
    return
end
if nargin==1
    parflag=false;
end

sol=odestruct(ocgtrajectory(ocgAsym),parflag);
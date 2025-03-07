function out=octype(ocgTrj)
% octype returns the type of the ocgradtrajectory, possible types are: 
% loc ... local (non-distributed) variables
% age ... age distributed
% space ... spatial distributed


deg=multiplicity(ocgTrj);
if deg>1
    out=cell(1,deg);
    for ii=1:deg
        out{ii}=ocgTrj(ii).type;
    end
else
    out=ocgTrj.type;
end

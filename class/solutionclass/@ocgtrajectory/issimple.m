function b=issimple(ocgTrj)
% ISSIMPLE checks if a generalized octrajectory can also be written as a
% simple octrajectory.
%
% OCGTRJ is simple if the number of ODEs is equal for all arcs.

b=false;
if isempty(ocgTrj)
    b=true;
    return
end
if all(ocgTrj.odenum==ocgTrj.odenum(1))
    b=true;
end
function ocTrj=octrajectory(ocgTrj)
%
% OCTRAJECTORY transforms a generalized octrajectory into an octrajectory if OCGTRJ is simple. 

if isempty(ocgTrj)
    ocTrj=octrajectory();
    return
end
if ~issimple(ocgTrj)
    ocmatmsg('Generalized ''octrajectory'' is not a simple ''octrajectory''.')
    ocTrj=octrajectory();
    return
end
    
ocTrj=ocgTrj.octrajectory;
ocTrj.y=[ocTrj.y{:}];

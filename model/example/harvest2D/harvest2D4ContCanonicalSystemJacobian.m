function [J,Jpar] = harvest2D4ContCanonicalSystemJacobian(t,dynVar,pararg,arcarg,switchid,switchcoord,trunct,internalpararg)
%
% the dynamics file for model harvest2D for the BVP
                                                                                                                                        
% Jacobian system

switch length(switchid)
    case 0
        J=trunct*harvest2DCanonicalSystemJacobian([],dynVar(1:4,:),pararg,arcarg(1));
        Jpar=zeros(4,1);
        %ADDCASE
    otherwise
        error([num2str(length(switchid)) ' %s not defined'],switchid)
end

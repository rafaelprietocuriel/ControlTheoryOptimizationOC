function out = harvest2D4ContCanonicalSystem(t,dynVar,pararg,arcarg)
%
% the dynamics file for model harvest2D for the BVP

% Dynamical system
% global variable

out(1:4,:)=harvest2DCanonicalSystem(t,dynVar(1:4,:),pararg,arcarg);

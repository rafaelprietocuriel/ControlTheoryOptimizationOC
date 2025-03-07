function out=dependentvar(gdynPrim)
%
% DEPENDENTVAR returns the point(s) of the dynprimitive object.
%
% P=DEPENDENTVAR(DYNPRIM) DYNPRIM is a member of the class
% 'dynprimitive', these can be interpreted either as a geometric object,
% point(s) in the phase space, or a time trajectory. The constant function
% in the case of an equilibrium, a periodic function in the case of a limit
% cycle. From the latter interpretation DEPENDENTVAR returns the (time)
% depending point(s) P of the object.


out=dependentvar(gdynPrim.ocgtrajectory);
function realdynPrim=real(dynPrim)
%
% REAL returns the real valued object.
%
% RDYNPRIM = REAL(DYNPRIM) returns a dynprimitive object RDYNPRIM, with only 
% the real parts, as well for the dynVar field as the linearization field,
% of the dynprimitive DYNPRIM.

realdynPrim=dynPrim;
realdynPrim.octrajectory=real(realdynPrim.octrajectory);
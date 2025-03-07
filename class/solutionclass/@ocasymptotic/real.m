function realocAsym=real(ocAsym)
%
% REAL returns the real valued counterpart.
%
% RDYNPRIM = REAL(DYNPRIM) returns a dynprimitive object RDYNPRIM, with only 
% the real parts, as well for the dynVar field as the linearization field,
% of the dynprimitive DYNPRIM.

realocAsym=ocAsym;
realocAsym.octrajectory=real(realocAsym.octrajectory);
realocAsym.limitset=real(realocAsym.limitset);
function out=arcargument(dynPrim)
%
% ARCARGUMENT returns the identfier for the underlying combination of
% active and inactive constraints.
%
% ARCARG=ARCARGUMENT(DYNPRIM) DYNPRIM is a member of the class
% 'dynprimitive'. The returned argument ARCARG is the identifier specified
% within an OCMat model for a specific combination of active (binding) and
% inactive (non-binding) constraints.

out=arcargument(dynPrim.octrajectory);
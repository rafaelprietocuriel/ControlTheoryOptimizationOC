function realdynPrim=remspuriousimag(dynPrim,opt)
%
% REMSPURIOUSIMAG removes spurious imaginary part.
%
% RDYNPRIM = REMSPURIOUSIMAG(DYNPRIM) returns a dynprimitive object
% RDYNPRIM, with only  the real parts, as well for the dynVar field as the
% linearization field, of the dynprimitive DYNPRIM.

if nargin==1
    opt=defaultocoptions;
end

if max(abs(imag(dynPrim)))>getocoptions(opt,'GENERAL','ImaginaryTolerance')
    realdynPrim=dynPrim;
    return
else
    realdynPrim=real(dynPrim);
end
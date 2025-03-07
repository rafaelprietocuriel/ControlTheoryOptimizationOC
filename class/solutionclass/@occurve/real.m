function realocCuv=real(ocCuv)
%
% IMAG returns an occurve object consisting of the imaginary parts.
%
% IOCCUV = IMAG(OCCUV) returns a dynprimitive object IOCCUV, with only 
% the imaginary parts, as well for the dynVar field as the linearization field,
% of the dynprimitive DYNPRIM.

realocCuv=ocCuv;
realocCuv.y=real(realocCuv.y);

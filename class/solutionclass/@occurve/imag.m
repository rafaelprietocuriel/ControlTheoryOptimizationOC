function varargout=imag(ocCuv)
%
% IMAG returns an occurve object consisting of the imaginary parts.
%
% IOCCUV = IMAG(OCCUV) returns a dynprimitive object IOCCUV, with only 
% the imaginary parts, as well for the dynVar field as the linearization field,
% of the dynprimitive DYNPRIM.

if nargout==0
    varargout{1}=imag(dependentvar(ocCuv));
else
    [varargout{1:nargout}]=imag(dependentvar(ocCuv));
end
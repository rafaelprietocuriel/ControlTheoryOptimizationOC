function varargout=imag(ppdePrim)
%
% IMAG returns a dynprimitive object consisting of the imaginary parts.
%
% IDYNPRIM = IMAG(DYNPRIM) returns a dynprimitive object RDYNPRIM, with only 
% the imaginary parts, as well for the dynVar field as the linearization field,
% of the dynprimitive DYNPRIM.

if nargout==0
    varargout{1}=imag(dependentvar(ppdePrim));
else
    [varargout{1:nargout}]=imag(dependentvar(ppdePrim));
end
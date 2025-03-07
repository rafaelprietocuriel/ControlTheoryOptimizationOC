function [b,dim]=issaddle(varargin)
%
% ISSADDLE true if limit set is of saddle type
%
% B = ISSADDLE(DYNPRIM) the input argument DYNPRIM is a dynprimitve object,
% 'equilibrium' or 'limitcycle'. B is 1 if DYNPRIM is of saddle type, for
% an equilibrium at least one eigenvalue has to be negative. For a limit
% cycle the absolute value at least of one eigenvalue has to be smaller
% than one.
%
% [B DIM] = ISSADDLE(DYNPRIM) DIM returns the number of of "stable"
% eigenvalues
%
% [B DIM] = ISSADDLE(DYNPRIM,OPT) with OPT the tolerance for a zero
% eigenvalue can be provided:
%   opt=setocoptions('OC','ZeroEigenValueTolerance',tol);
%
% [B DIM] = ISSADDLE(DYNPRIM1,...,DYNPRIMN) B is a vector of 0 and 1
% depending if the dynprimitve objects DYNPRIM1,...,DYNPRIMN are of saddle
% type. DIM is a vector containing the number of "stable" eigenvlaues for
% each of the dynprimitve objects.
%
% [B DIM] = ISSADDLE(DYNPRIM1,...,DYNPRIMN,OPT)

[numseigval,numueigval]=characteristics(varargin{:});
dim=[numseigval{:}];
b=dim & [numueigval{:}]>0;

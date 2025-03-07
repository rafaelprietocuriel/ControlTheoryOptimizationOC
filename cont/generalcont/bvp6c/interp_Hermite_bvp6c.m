function [Sx, Spx] = interp_Hermite_bvp6c(w,h,y,yp,yp_ip025,yp_ip05,yp_ip075)
%INTERP_HERMITE  use the 6th order Hermite Interpolant presented by Cash
%and Wright to find y and y' at abscissas for Quadrature.
N=size(y,2);
diagscal = spdiags(h',0,N-1,N-1);

Sx=   A66(w)*y(:,2:N) + A66(1-w)*y(:,1:N-1)   + ...
    ( B66(w)*yp(:,2:N) - B66(1-w)*yp(:,1:N-1) + ...
    C66(w)*(yp_ip075-yp_ip025) + D66(w)*yp_ip05 )*diagscal;

if nargout >1
    diagscal = spdiags(1./h',0,N-1,N-1);
    Spx=( Ap66(w)*y(:,2:N)  - Ap66(1-w)*y(:,1:N-1) )*diagscal + ...
        ( Bp66(w)*yp(:,2:N) + Bp66(1-w)*yp(:,1:N-1) + ...
        Cp66(w)*(yp_ip075-yp_ip025) + Dp66(w)*yp_ip05 );
end

%---------------------------------------------------------------------------
function coeff = A66(w)
coeff=w.^2.*polyval([-24 60 -50 15],w);     % w^2*(15-50*w+60*w^2-24*w^3);
function coeff = B66(w)
coeff=w.^2.*polyval([12 -26 19 -5]/3,w);    % w^2/3*(w-1)*(12*w^2-14*w+5);
function coeff = C66(w)
coeff=w.^2.*polyval([-8 16 -8]/3,w);        % -w^2*8/3*(1-w)^2;
function coeff = D66(w)
coeff=w.^2.*polyval([16 -40 32 -8],w);      % w^2*8*(1-w)^2*(2*w-1);

function coeff = Ap66(w)
coeff=w.*polyval([-120 240 -150 30],w);     %w*(30-150*w+240*w^2-120*w^3);
function coeff = Bp66(w)
coeff=w.*polyval([20 -104/3 19 -10/3],w);   %w*(w*(20*w^2+19)-(104*w^2+10)/3);
function coeff = Cp66(w)
coeff=-16/3*w.*polyval([2 -3 1],w);         %-16/3*w*(1-3*w+2*w^2);
function coeff = Dp66(w)
coeff=w.*polyval([80 -160 96 -16],w);       %w*(80*w^3-160*w^2+96*w-16);

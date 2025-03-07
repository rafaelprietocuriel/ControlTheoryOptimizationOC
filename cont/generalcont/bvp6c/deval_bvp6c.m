function [Sxint,Spxint] = deval_bvp6c(t,y,xint,yp,ypmid)
%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2004 The MathWorks, Inc.
%   BVP6C Modification
%    Nick Hale  Imperial College London
%    $Date: 12/06/2006 $

% Determine the dominant data type
dataType = superiorfloat(t,xint);

Spxint_requested = (nargout > 1);   % Evaluate the first derivative?

n = size(y,1);
Nxint = length(xint);
Sxint = zeros(n,Nxint,dataType);
if Spxint_requested
  Spxint = zeros(n,Nxint,dataType);
end

% Make tint a row vector and if necessary, 
% sort it to match the order of t.
tint = xint(:).';  
tdir = sign(t(end) - t(1));
had2sort = any(tdir*diff(tint) < 0);
if had2sort
  [tint,tint_order] = sort(tdir*tint);
  tint = tdir*tint;
end  

% Using the sorted version of tint, test for illegal values.
if (tdir*(tint(1) - t(1)) < 0) || (tdir*(tint(end) - t(end)) > 0)
  ocmaterror('MATLAB:deval:SolOutsideInterval',...
        ['Attempting to evaluate the solution outside the interval\n'...
         '[%e, %e] where it is defined.\n'],t(1),t(end));
end
interpfcn = @ntrp6h;

evaluated = 0;
bottom = 1;
while evaluated < Nxint
  
  % Find right-open subinterval [t(bottom), t(bottom+1))
  % containing the next entry of tint. 
  Index = find( tdir*(tint(evaluated+1) - t(bottom:end)) >= 0 );
  bottom = bottom - 1 + Index(end);

  % Is it [t(end), t(end)]?
  at_tend = (t(bottom) == t(end));
      
  if at_tend
    % Use (t(bottom-1) t(bottom)] to interpolate y(t(end)) and yp(t(end)).
    index = find(tint(evaluated+1:end) == t(bottom));
    bottom = bottom - 1;    
  else    
    % Interpolate inside [t(bottom) t(bottom+1)).
    index = find( tdir*(tint(evaluated+1:end) - t(bottom+1)) < 0 );  
  end

  interpdata = { yp(:,bottom) yp(:,bottom+1) ...
        ypmid(:,bottom,:)};               
  
  % Evaluate the interpolant at all points from [t(bottom), t(bottom+1)).
  if Spxint_requested
    [yint,ypint] = feval(interpfcn,tint(evaluated+index),t(bottom),y(:,bottom),...
                         t(bottom+1),y(:,bottom+1),interpdata{:});    
  else  
    yint = feval(interpfcn,tint(evaluated+index),t(bottom),y(:,bottom),...
                 t(bottom+1),y(:,bottom+1),interpdata{:});    
  end

  if at_tend
    bottom = bottom + 1;
  end
  
  % Purify the solution at t(bottom).
  index1 = find(tint(evaluated+index) == t(bottom));
  if ~isempty(index1)    
    yint(:,index1) = repmat(y(:,bottom),1,length(index1)); 
  end
  
  % Accumulate the output.
  Sxint(:,evaluated+index) = yint(:,:);  
  if Spxint_requested
    Spxint(:,evaluated+index) = ypint(:,:);
  end  
  evaluated = evaluated + length(index);  
end

% For multipoint BVPs, check if solution requested at interface points
multipointBVP = (any(diff(t)==0));
if multipointBVP
  idxDiscontPoints = find(diff(t)==0);
  for i = idxDiscontPoints
    % Check whether solution requested at the interface
    idxIntPoints = find(tint == t(i));
    n = length(idxIntPoints);
    if n > 0
      % Average the solution across the interface
%       nl = '\n         ';
%       warning('MATLAB:deval:NonuniqueSolution',...
%               ['At the interface XC = %g the solution might not be continuous.'...
%                nl 'DEVAL returns the average of the limits from the left and'...
%                nl 'from the right of the interface. To approximate the limit'... 
%                nl 'values, call DEVAL for XC-EPS(XC) or XC+EPS(XC).\n'],t(i));
      
      Sxint(:,idxIntPoints) = repmat((y(:,i)+y(:,i+1))/2,1,n);
      if Spxint_requested
        Spxint(:,idxIntPoints) = repmat((yp(:,i)+yp(:,i+1))/2,1,n);
      end  
    end    
  end  
end
  
if had2sort     % Restore the order of tint in the output.
  Sxint(:,tint_order) = Sxint;  
  if Spxint_requested
    Spxint(:,tint_order) = Spxint;  
  end  
end

function [yint,ypint] = ntrp6h(tint,t,y,tnew,ynew,yp,ypnew,Fmid)
%NTRP6H  Interpolation helper function for BVP6C.
%   YINT = NTRP6H(TINT,T,Y,TNEW,YNEW,YP,YPNEW,FMID) evaluates the
%   Cash-Moore, Cash-Singhal6 based linear interpolant at time TINT.
%   TINT may be a scalar or a row vector.   
%   [YINT,YPINT] = NTRP6H(TINT,T,Y,TNEW,YNEW,YP,YPNEW,FMID) returns
%   also the derivative of the interpolating polynomial. 
%   
%   See also BVP6C, DEVAL, NTRP6C

%    Nick Hale  Imperial College London
%    $Date: 12/06/2006 $

h = tnew - t;
w = (tint - t)/h;
for i=1:length(tint)
    yint(:,i) = A66(w(i))*ynew + A66(1-w(i))*y   + ...
            ( B66(w(i))*ypnew - B66(1-w(i))*yp + ...
              C66(w(i))*(Fmid(:,:,3)-Fmid(:,:,1)) + D66(w(i))*Fmid(:,:,2) )*h;         
    if nargout > 1
        ypint(:,i) =( Ap66(w(i))*ynew  - Ap66(1-w(i))*y )/h + ...
               ( Bp66(w(i))*ypnew + Bp66(1-w(i))*yp + ...
                 Cp66(w(i))*(Fmid(:,:,3)-Fmid(:,:,1)) + Dp66(w(i))*Fmid(:,:,2) );  
    end
end    

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


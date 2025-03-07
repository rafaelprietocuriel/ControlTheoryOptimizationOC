function [Sxint,Spxint] = deval_bvp5c(t,y,xint,yp,ymid)
%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2004 The MathWorks, Inc.
%   BVP6C Modification
%    Nick Hale  Imperial College London
%    $Date: 12/06/2006 $

% Determine the dominant data type
dataType = superiorfloat(t,xint);

Spxint_requested = (nargout > 1);   % Evaluate the first derivative?

sol.x=t;
sol.y=y;
sol.idata.yp=yp;
sol.idata.ymid=ymid;
[n,N]=size(y);
idx = 1:n;
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
interpfcn = @ntrp4h;  

evaluated = 0;
bottom = 1;
while evaluated < Nxint

    % Find right-open subinterval [t(bottom), t(bottom+1))
    % containing the next entry of tint.
    Index = find( tdir*(tint(evaluated+1) - t(bottom:end)) >= 0 );
    bottom = bottom - 1 + Index(end);

    % Is it [t(end), t(end)]?
    at_tend = (t(bottom) == t(end));

    % Return solution already available at t(bottom)
    index1 = find(tint(evaluated+1:end) == t(bottom));

    % Interpolate solution inside (t(bottom), t(bottom+1))
    if at_tend
        index2 = [];
    else
        index2 = find( (tdir*(tint(evaluated+1:end) - t(bottom)) > 0) & ...
            (tdir*(tint(evaluated+1:end) - t(bottom+1)) < 0) );
    end

    % Return the (adjusted) solution at t(bottom)
    if ~isempty(index1)
        if at_tend
            yint1 = y(:,end);
            if Spxint_requested
                % Extrapolate derivative from [t(bottom-1),t(bottom))
                interpdata = extract_idata(sol,t,bottom-1);
                [ignore,ypint1] = interpfcn(t(bottom),t(bottom-1),y(:,bottom-1),...
                    t(bottom),y(:,bottom),interpdata{:});
            end

        elseif (bottom > 2) && (t(bottom) == t(bottom-1)) % Interface point
            % Average the solution (and its derivative) across the interface.
            yLeft  = y(:,bottom-1);
            yRight = y(:,bottom);
            yint1 = (yLeft + yRight)/2;
            if Spxint_requested
                % Get the 'left' derivative by extrapolating in [t(bottom-2), t(bottom-1)).
                interpdata =  extract_idata(sol,t,bottom-2);
                [ignore,ypLeft] = interpfcn(t(bottom-1),t(bottom-2),y(:,bottom-2),...
                    t(bottom-1),y(:,bottom-1),interpdata{:});
                % Get the 'right' derivative by interpolating in [t(bottom), t(bottom+1)).
                interpdata =  extract_idata(sol,t,bottom);
                [ignore,ypRight] = interpfcn(t(bottom),t(bottom),y(:,bottom),...
                    t(bottom+1),y(:,bottom+1),interpdata{:});
                ypint1 = (ypLeft + ypRight)/2;
            end
%             nl = '\n         ';
%             warning('MATLAB:deval:NonuniqueSolution',...
%                 ['At the interface XC = %g the solution might not be continuous.'...
%                 nl 'DEVAL returns the average of the limits from the left and'...
%                 nl 'from the right of the interface. To approximate the limit'...
%                 nl 'values, call DEVAL for XC-EPS(XC) or XC+EPS(XC).\n'],t(bottom));
% 
        else
            % 'Regular' mesh point
            yint1 = y(:,bottom);
            if Spxint_requested
                % Interpolate derivative from [t(bottom),t(bottom+1))
                interpdata = extract_idata(sol,t,bottom);
                [ignore,ypint1] = interpfcn(t(bottom),t(bottom),y(:,bottom),...
                    t(bottom+1),y(:,bottom+1),interpdata{:});
            end
        end

        % Accumulate the output.
        Sxint(:,evaluated+index1) = yint1(idx,ones(1,numel(index1)));
        if Spxint_requested
            Spxint(:,evaluated+index1) = ypint1(idx,ones(1,numel(index1)));
        end
    end

    % Interpolate solution inside (t(bottom), t(bottom+1)).
    if ~isempty(index2)
        % Get 'bvp4c'-dependent interpolation data for [t(bottom), t(bottom+1)).
        interpdata = extract_idata(sol,t,bottom);

        % Evaluate the interpolant at all points from (t(bottom), t(bottom+1)).
        if Spxint_requested
            [yint2,ypint2] = interpfcn(tint(evaluated+index2),t(bottom),y(:,bottom),...
                t(bottom+1),y(:,bottom+1),interpdata{:});
        else
            yint2 = interpfcn(tint(evaluated+index2),t(bottom),y(:,bottom),...
                t(bottom+1),y(:,bottom+1),interpdata{:});
        end

        % Accumulate the output.
        Sxint(:,evaluated+index2) = yint2(idx,:);
        if Spxint_requested
            Spxint(:,evaluated+index2) = ypint2(idx,:);
        end
    end

    evaluated = evaluated + length(index1) + length(index2);
end

if had2sort     % Restore the order of tint in the output.
    Sxint(:,tint_order) = Sxint;
    if Spxint_requested
        Spxint(:,tint_order) = Spxint;
    end
end

% --------------------------------------------------------------------------

function interpdata = extract_idata(sol,t,tidx)
% Data for interpolation in [t(tidx), t(tidx+1))

  interpdata = { sol.idata.ymid(:,tidx), ...
                 sol.idata.yp(:,tidx), ...
                 sol.idata.yp(:,tidx+1) };    

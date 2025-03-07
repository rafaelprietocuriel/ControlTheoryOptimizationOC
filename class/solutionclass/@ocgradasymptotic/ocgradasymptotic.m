function ocgAsym=ocgradasymptotic(varargin)
%
% OCGRADASYMPTOTIC ocgradasymptotic constructor
%
% OCGRADASYMPTOTIC(ODESTRUCT) creates an ocgradasymptotic object from an ode solution
% structure ODESTRUCT.
%
switch nargin
    case 0
        ocgAsym.limitset=limitset(dynprimitive());
        ocgAsym=class(ocgAsym,'ocgradasymptotic',ocgradtrajectory());
        
    case 1
        if isocgradasymptotic(varargin{1})
            if isempty(varargin{1})
                ocgAsym=ocgradasymptotic();
            else
                ocgAsym=varargin{1};
            end
        elseif isempty(varargin{1})
            ocgAsym=ocgradasymptotic();
        end
    case 2
        if isocgradtrajectory(varargin{1}) || icocgradlimitset(varargin{2})
            ocgAsym.limitset=varargin{2};
            ocgAsym=class(ocgAsym,'ocgradasymptotic',varargin{1});
        elseif isocgradtrajectory(varargin{2}) || icocgradlimitset(varargin{1})
            ocgAsym.limitset=varargin{1};
            ocgAsym=class(ocgAsym,'ocgradasymptotic',varargin{2});
        else
            ocgAsym=ocgradasymptotic();
        end
end

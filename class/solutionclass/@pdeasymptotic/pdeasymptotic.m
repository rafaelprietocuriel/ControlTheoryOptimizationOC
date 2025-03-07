function pdeAsym=pdeasymptotic(varargin)
%
% PDEASYMPTOTIC pdeasymptotic constructor
%
switch nargin
    case 0
        pdeAsym.pathtype='';
        %pdeAsym.pdeprimitive=pdeprimitive();
        pdeAsym=class(pdeAsym,'pdeasymptotic',pdetrajectory(),pdeprimitive());
        
    case 1
        if ispdeasymptotic(varargin{1})
            pdeAsym=pdeasymptotic();
        elseif isempty(varargin{1})
            pdeAsym=pdeasymptotic();
        elseif ispdeprimitive(varargin{1})
            % an equilibrium (limit cycle) as an asymptotic solution
            pdeAsym.pathtype='';
            pdeAsym=class(pdeAsym,'pdeasymptotic',pdetrajectory(varargin{1}),varargin{1});
        end
    case 2
        if ispdeprimitive(varargin{2})
            pdePrim=varargin{2};
            if ispdetrajectory(varargin{1}) || ispdeasymptotic(varargin{1})
                pdeTrj=pdetrajectory(varargin{1});
            elseif isstruct(varargin{1})
                pdeTrj=pdetrajectory(varargin{1});
            end
            % an equilibrium (limit cycle) as an asymptotic solution
            pdeAsym.pathtype='';
            pdeAsym=class(pdeAsym,'pdeasymptotic',pdeTrj,pdePrim);
        else
            ocmatmsg('Second argument is not a ''pdeprimitive''.')
            pdeAsym=pdeasymptotic();
        end
    case 3
        if ispdetrajectory(varargin{1}) && ispdeprimitive(varargin{2}) && ischar(varargin{3})
            % an equilibrium (limit cycle) as an asymptotic solution
            pdeAsym.pathtype=varargin{3};
            pdeAsym=class(pdeAsym,'pdeasymptotic',varargin{1},varargin{2});
        else
            ocmatmsg('Initialization failed.')
            pdeAsym=pdeasymptotic();
        end
        
end
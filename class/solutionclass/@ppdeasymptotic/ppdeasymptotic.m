function ppdeAsym=ppdeasymptotic(varargin)
%
% PPDEASYMPTOTIC ppdeasymptotic constructor
%
switch nargin
    case 0
        ppdeAsym.limitset=ppdeprimitive();
        ppdeAsym=class(ppdeAsym,'ppdeasymptotic',ppdetrajectory());
        
    case 1
        if isppdeasymptotic(varargin{1})
            if isempty(varargin{1})
                ppdeAsym=ppdeasymptotic();
            else
                ppdeAsym=varargin{1};
            end
        elseif isempty(varargin{1})
            ppdeAsym=ppdeasymptotic();
        elseif isppdeprimitive(varargin{1})
            % an equilibrium (limit cycle) as an asymptotic solution
            ppdeAsym.limitset=varargin{1};
            ppdeAsym=class(ppdeAsym,'ppdeasymptotic',ppdetrajectory(ppdeAsym));
        end
    case 2
        try
            % first argument ppdetrajectory
            % second argument ppdeprimitive
            ppdePrim=varargin{2};
            ppdeAsym.limitset=ppdePrim;
            ppdeAsym=class(ppdeAsym,'ppdeasymptotic',varargin{1});
            return
        end
        try
            % first argument ppdeprimitive
            % second argument ppdetrajectory
            ppdePrim=ppdeprimitive(varargin{1});
            ppdeAsym.limitset=ppdePrim;
            ppdeAsym=class(ppdeAsym,'ppdeasymptotic',ppdetrajectory(varargin{2}));
        catch
            lasterr
        end
end
function ocAsym=docasymptotic(varargin)
%
% OCASYMPTOTIC ocasymptotic constructor
%
switch nargin
    case 0
        ocAsym.limitset=limitset(mapprimitive());
        ocAsym=class(ocAsym,'docasymptotic',doctrajectory(),mapprimitive());
        
    case 1
        if isdocasymptotic(varargin{1})
            if isempty(varargin{1})
                ocAsym=docasymptotic();
            else
                ocAsym=varargin{1};
            end
        elseif isempty(varargin{1})
            ocAsym=docasymptotic();
        elseif ismapprimitive(varargin{1})
            % an equilibrium (limit cycle) as an asymptotic solution
            ocAsym.limitset=limitset(varargin{1});
            ocAsym=class(ocAsym,'docasymptotic',octrajectory(ocAsym));
        end
    case 2
        try
            % first argument octrajectory
            % second argument mapprimitive
            dynPrim=mapprimitive(varargin{2});
            ocAsym.limitset=limitset(dynPrim);
            ocAsym=class(ocAsym,'docasymptotic',doctrajectory(varargin{1}),dynPrim);
            return
        end
        try
            % first argument mapprimitive
            % second argument doctrajectory
            dynPrim=mapprimitive(varargin{1});
            ocAsym.limitset=limitset(dynPrim);
            ocAsym=class(ocAsym,'docasymptotic',doctrajectory(varargin{2}),dynPrim);
        catch
            lasterr
        end
end
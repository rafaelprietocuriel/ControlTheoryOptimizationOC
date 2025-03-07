function sol=generatesolstruct(ocAsym,solvername,varargin)

sol=generatesolstruct(ocAsym.octrajectory,solvername,varargin{:});
sol.timehorizon=inf;
function ocAsym=redefinearc(ocAsym,newposition,arcid,varargin)

ocAsym=ocasymptotic(redefinearc(octrajectory(ocAsym),newposition,arcid,varargin{:}),ocAsym.dynprimitive);
function ocgAsym=redefinearc(ocgAsym,newposition,arcid,varargin)

ocgAsym=ocgasymptotic(redefinearc(ocgtrajectory(ocgAsym),newposition,arcid,varargin{:}),limitset(ocgAsym));
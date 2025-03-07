function ocTrj=redefinearc(ocTrj,newposition,arcid,varargin)
pathtpe=pathtype(ocTrj);
ocTrj=octrajectory(redefinearc(struct(ocTrj),newposition,arcid,varargin{:}));

ocTrj.solverinfo.pathtype=pathtpe;
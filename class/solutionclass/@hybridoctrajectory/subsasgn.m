function ocTrj=subsasgn(ocTrj,index,rhs)
%
%
try
    ocTrj=hybridoctrajectory(subsasgn(struct(ocTrj),index,rhs));
catch
    lasterr;
end
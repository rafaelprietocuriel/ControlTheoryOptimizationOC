function ocTrj=subsasgn(ocTrj,index,rhs)
%
%
try
    ocTrj=pdetrajectory(subsasgn(struct(ocTrj),index,rhs));
catch
    lasterr;
end
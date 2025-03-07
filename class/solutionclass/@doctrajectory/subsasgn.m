function docTrj=subsasgn(docTrj,index,rhs)
%
%
try
    docTrj=doctrajectory(subsasgn(struct(docTrj),index,rhs));
catch
    lasterr;
end
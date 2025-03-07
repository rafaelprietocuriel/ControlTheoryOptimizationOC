function docAsym=subsasgn(docAsym,index,rhs)
%
%
try
    o=subsasgn(struct(docAsym),index,struct(rhs));
catch
    o=subsasgn(struct(docAsym),index,rhs);
end
docAsym=docasymptotic([o.doctrajectory],[o.limitset]);

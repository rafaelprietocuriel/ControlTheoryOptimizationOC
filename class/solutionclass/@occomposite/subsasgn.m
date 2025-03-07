function ocComp=subsasgn(ocComp,index,rhs)
%
%
try
    ocComp=occomposite(subsasgn(struct(ocComp),index,rhs));
catch
    lasterr;
end
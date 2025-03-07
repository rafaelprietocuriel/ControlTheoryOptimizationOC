function ocCuv=subsasgn(ocCuv,index,rhs)
%
%
try
    ocCuv=subsasgn(struct(ocCuv),index,rhs);
catch
    lasterr;
end
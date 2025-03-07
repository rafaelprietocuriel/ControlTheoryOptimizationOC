function mapPrim=subsasgn(mapPrim,index,rhs)
%
%
try
    if strcmp(index.subs,'period')
        mapPrim=mapprimitive(subsasgn(struct(mapPrim),index,rhs));
    else
        mapPrim.doctrajectory=subsasgn(mapPrim.doctrajectory,index,rhs);
    end
catch
    lasterr;
end
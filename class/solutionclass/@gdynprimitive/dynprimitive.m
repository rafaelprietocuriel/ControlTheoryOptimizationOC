function dynPrim=dynprimitive(gdynPrim)
%
%

if issimple(gdynPrim.ocgtrajectory)
    dynPrim.octrajectory=octrajectory(gdynPrim.ocgtrajectory);
    dynPrim.period=gdynPrim.period;
    dynPrim=dynprimitive(dynPrim);
else
    dynPrim=dynprimitive();
end

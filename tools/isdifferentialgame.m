function b=isdifferentialgame(ocObj)

b=~isempty(regexp(class(ocObj),'differentialgame\>','ONCE'));
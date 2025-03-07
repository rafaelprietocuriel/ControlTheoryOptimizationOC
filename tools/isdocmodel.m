function b=isdocmodel(ocObj)

b=~isempty(regexp(class(ocObj),'docmodel\>','ONCE'));
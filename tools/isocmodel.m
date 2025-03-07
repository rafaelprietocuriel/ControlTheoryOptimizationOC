function b=isocmodel(ocObj)

b=~isempty(regexp(class(ocObj),'ocmodel\>','ONCE'));
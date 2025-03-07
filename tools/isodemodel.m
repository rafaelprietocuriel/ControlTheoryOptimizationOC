function b=isodemodel(ocObj)

b=~isempty(regexp(class(ocObj),'\<odemodel\>','ONCE'));
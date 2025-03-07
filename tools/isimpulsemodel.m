function b=isimpulsemodel(ocObj)

b=~isempty(regexp(class(ocObj),'\<impulsemodel\>','ONCE'));
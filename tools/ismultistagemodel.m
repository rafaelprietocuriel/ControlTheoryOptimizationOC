function b=ismultistagemodel(ocObj)

b=~isempty(regexp(class(ocObj),'multistagemodel\>','ONCE'));
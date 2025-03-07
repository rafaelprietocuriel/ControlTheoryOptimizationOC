function b=ismultimodel(mmObj)

b=~isempty(regexp(class(mmObj),'multimodel\>','ONCE'));
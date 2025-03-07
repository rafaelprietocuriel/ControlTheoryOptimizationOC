function b=isppdemodel(ppdeObj)

b=~isempty(regexp(class(ppdeObj),'ppdemodel\>','ONCE'));
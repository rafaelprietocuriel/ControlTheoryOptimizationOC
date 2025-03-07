function wsn=wsname(ppdeocObj,number)
%
% workspace name depending on oc model

wsn='';
if isempty(ppdeocObj)
    return
end
if nargin==1
    number='';
else
    if ~ischar(number)
        number=num2str(number);
    end
end
wsn=fullocmatfile(userdatafolder(ppdeocObj),[modelname(ppdeocObj) 'WorkSpace' number '.mat']);
function wsn=wsname(mmObj,number)
%
% workspace name depending on oc model

wsn='';
if isempty(mmObj)
    return
end
if nargin==1
    number='';
else
    if ~ischar(number)
        number=num2str(number);
    end
end
wsn=fullocmatfile(userdatafolder(mmObj),[modelname(mmObj) 'WorkSpace' number '.mat']);
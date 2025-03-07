function wsn=wsname(ocObj,number)
%
% workspace name depending on oc model

wsn='';
if isempty(ocObj)
    return
end
if nargin==1
    number='';
else
    if ~ischar(number)
        number=num2str(number);
    end
end
wsn=fullocmatfile(userdatafolder(ocObj),[modelname(ocObj) 'WorkSpace' number '.mat']);
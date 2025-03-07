function wsn=wsname(dgObj,number)
%
% workspace name depending on oc model

wsn='';
if isempty(dgObj)
    return
end
if nargin==1
    number='';
else
    if ~ischar(number)
        number=num2str(number);
    end
end
wsn=fullocmatfile(userdatafolder(dgObj),[modelname(dgObj) 'WorkSpace' number '.mat']);
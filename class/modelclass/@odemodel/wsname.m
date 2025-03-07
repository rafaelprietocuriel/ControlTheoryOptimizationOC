function wsn=wsname(odeObj,number)
%
% workspace name depending on oc model

wsn='';
if isempty(odeObj)
    return
end
if nargin==1
    number='';
else
    if ~ischar(number)
        number=num2str(number);
    end
end
wsn=fullocmatfile(userdatafolder(odeObj),[modelname(odeObj) 'WorkSpace' number '.mat']);
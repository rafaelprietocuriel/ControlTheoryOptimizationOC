function showparametervalue(dgObj,varargin)
%
% SHOWPARAMETERVALUE displays the parameter value(s) of the oc model
%
% SHOWPARAMETERVALUE(OCOBJ)
%
% SHOWPARAMETERVALUE(OCOBJ,FORMAT) displays the parameter values, where
% FORMAT can either be a number, denoting the maximum precision for the
% parameter values or a string denoting the format. (See fprintf for format
% string details.) 

if isempty(dgObj)
    ocmatmsg('Oc model is empty.\n')
    return
end
num2strarg=[];

if nargin==2
    num2strarg=varargin{1};
end

if isempty(num2strarg)
    num2strarg='%3.4g';
end
[parval,parvar]=parametervalue(dgObj);

for ii=1:parameternum(dgObj)
    disp([parvar{ii} ' : ' num2str(parval(ii),num2strarg)]);
end
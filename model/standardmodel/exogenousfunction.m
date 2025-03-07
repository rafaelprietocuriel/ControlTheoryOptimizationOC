function b=exogenousfunction(ocStruct)
%
%
b=[];
if isempty(ocStruct)
    return
end
out=retrievemodelinformation(ocStruct,'exogenousfunctionnum');
b=out.value~=0;
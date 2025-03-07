function b=exogenousdynamics(ocStruct)
%
%
b=[];
if isempty(ocStruct)
    return
end
out=retrievemodelinformation(ocStruct,'exogenousdynamicsnum');
b=out.value~=0;
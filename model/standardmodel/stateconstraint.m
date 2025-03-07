function b=stateconstraint(ocStruct)
%
%
b=[];
if isempty(ocStruct)
    return
end
out=retrievemodelinformation(ocStruct,'inequalitystateconstraintnum');
b=out.value~=0;
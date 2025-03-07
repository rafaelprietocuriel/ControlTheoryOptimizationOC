function b=nonsmoothfunction(ocStruct)
%
%
b=[];
if isempty(ocStruct)
    return
end
out=retrievemodelinformation(ocStruct,'nonsmoothfunctionnum');
b=out.value&&~isempty(out.value);
function [b,idx]=mixedconstraint(ocObj)
b=0;
idx=[];
if isempty(ocObj)
    return
end

cstr=constraint(ocObj);
x=state(ocObj);

idx=1:length(cstr);
for ii=1:length(cstr)
    idx(ii)=0;
    jj=0;
    while jj<length(x)
        jj=jj+1;
        if diff(cstr(ii),x{jj})~=0
            idx(ii)=ii;
            break
        end
    end
end
idx(idx==0)=[];
if any(idx)
    b=1;
end
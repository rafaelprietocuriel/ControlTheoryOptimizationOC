function out=algebraicequation(ocObj,arcid)
%
out=[];
if nargin==1
    arcid=0;
end
if isempty(ocObj)
    return
end

algebraicequation=retrievemodelinformation(ocObj,'algebraicequationimplicit',num2str(arcid),getsymkernel());
out=sym(cell2vectorstring(algebraicequation.value)).';

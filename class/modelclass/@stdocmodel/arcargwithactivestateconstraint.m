function arcid=arcargwithactivestateconstraint(ocObj)
%
% ARCIDENTIFIER arcidentifiers are characters of numbers, consecutively
% numerated.

arcid=[];
if isempty(ocObj)
    return
end
arcarg=arcargument(ocObj);
counter=0;
for ii=1:length(arcarg)
    info=retrievemodelinformation(ocObj.Model,'nonzerolmsc',num2str(arcarg(ii)));
    if ~isempty(info.value)
        counter=counter+1;
        arcid(counter)=arcarg(ii);
    end
end

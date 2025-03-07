function n=numberofparts(mmMultiPath)
%
%
if isempty(mmMultiPath)
    n=0;
else
    n=mmMultiPath.parts;
end
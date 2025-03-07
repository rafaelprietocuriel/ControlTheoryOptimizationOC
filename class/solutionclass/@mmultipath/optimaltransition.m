function n=optimaltransition(mmMultiPath)
%
%
if isempty(mmMultiPath)
    n=0;
else
    n=mmMultiPath.optimaltransition;
end
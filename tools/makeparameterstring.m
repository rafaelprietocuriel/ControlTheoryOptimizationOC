function parchar=makeparameterstring(parcell)

if ~iscell(parcell)
    parchar='';
    return
end

parchar='';
for ii=1:length(parcell)
    parchar=[parchar ',' parcell{ii}];
end
parchar(1)=[];


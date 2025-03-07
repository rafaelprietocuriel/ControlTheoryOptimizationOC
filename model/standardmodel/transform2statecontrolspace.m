function b=transform2statecontrolspace(ocStruct)
%
%
if isempty(ocStruct)
    b=false;
    return
end
modelclass=modeltype(ocStruct);
switch modelclass
    case 'standardmodel'
        arcidentifier=retrievemodelinformation(ocStruct,'arcidentifier');
        b=repmat(false,1,numel(arcidentifier.value));
        for ii=1:numel(arcidentifier.value)
            transformflag=retrievemodelinformation(ocStruct,'transform2statecontrolspace',arcidentifier.value{ii});
            b(ii)=transformflag.value==1;
        end
    case 'multistagemodel'
        b=false;
end
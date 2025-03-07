function extremal=generategradextremal(ocgMTrj)
% generate extremal for gradient initialization
extremal=[];
if isempty(ocgMTrj)
    return
end

ocgTrjC=ocgradmultipath2cell(ocgMTrj);
extremal=generategradextremal(ocgTrjC{1});
for ii=2:length(ocgTrjC)
    extremal(ii)=generategradextremal(ocgTrjC{ii});
end

function extremal=generatestaticextremal(ocgTrj,varargin)
% generate extremal for gradient initialization

extremal=[];
if isempty(ocgTrj)
    return
end
var=variable(ocgTrj);
fn=fieldnames(var);
for ii=1:length(fn)
    extremal.(fn{ii})=var.(fn{ii});
end
extremal.objectivevalue=objectivevalue(ocgTrj);
extremal.coeff=extremal.y(:);
extremal.dv=[];
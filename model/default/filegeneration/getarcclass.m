function arcid=getarcclass(ocStruct,class)
arcid='';

if isempty(ocStruct)
    return
end
switch class
    case 'arc'
        for ii=1:ocStruct.arc.num
            arcid{ii}=arcidentifier2field(ocStruct.arc.identifier{ii});
        end
    case {'algebraicequation'}
        % remove identifers where no algebraic equation exist
        arcid=fieldnames(ocStruct.foc.adjointsystem.algebraicequation);
    case {'daequation'}
        arcid=fieldnames(ocStruct.foc.adjointsystem.algebraicequation);
    case {'canonicalsystem'}
        arcid=fieldnames(ocStruct.foc.adjointsystem.dynamics);
    case {'control','controlvalue','canonicalsystemjacobian'}
        arcid=fieldnames(ocStruct.foc.value.control);
    case {'impulsecontrol','impulsecontrolvalue'}
        arcid=fieldnames(ocStruct.foc.value.icontrol);
    case {'lagrangemultcc','lagrangemultiplier'}
        if isfield(ocStruct.foc.value,'lagrangemultcc')
            arcid=fieldnames(ocStruct.foc.value.lagrangemultcc);
        elseif isfield(ocStruct.foc.value,'lagrangemultsc')
            arcid=getarcclass(ocStruct,'lagrangemultsc');
        end
    case {'lagrangemultsc'}
        arcid=fieldnames(ocStruct.foc.value.lagrangemultsc);
    case 'controlnum'
        arcid=1:ocStruct.variable.control.num;
end
if length(arcid)>9
    arcvalue=zeros(1,length(arcid));
    for ii=1:length(arcid)
        arcvalue(ii)=str2double(arcid{ii}(4:end));
    end
    [dum,idx]=sort(arcvalue);
    arcid=arcid(idx);
end

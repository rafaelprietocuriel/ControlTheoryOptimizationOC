function arcpos=getarcpositioncol(tcol)
% returns the position of the arcs for the collocation mesh for each stage

persistent stagenum N arcposition

if isempty(stagenum) || length(tcol)~=stagenum
    stagenum=length(tcol);
    for ii=1:stagenum
        N(ii)=length(tcol{ii});
        idx=find(diff(tcol{ii})==0);
        arcposition{ii}=[1 idx+1;idx N(ii)];
    end
    arcpos=arcposition;
    return
end

for ii=1:stagenum
    if N(ii)~=length(tcol{ii})
        N(ii)=length(tcol{ii});
        idx=find(diff(tcol{ii})==0);
        arcposition{ii}=[1 idx+1;idx N(ii)];
    end
end

arcpos=arcposition;

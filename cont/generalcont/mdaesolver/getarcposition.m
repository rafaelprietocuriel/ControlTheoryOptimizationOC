function arcpos=getarcposition(t)
% returns the position of the arcs for each stage
persistent stagenum N arcposition

if isempty(stagenum) || length(t)~=stagenum
    stagenum=length(t);
    for ii=1:stagenum
        N(ii)=length(t{ii});
        idx=find(diff(t{ii})==0);
        arcposition{ii}=[1 idx+1;idx N(ii)];
    end
    arcpos=arcposition;
    return
end

for ii=1:stagenum
    if N(ii)~=length(t{ii})
        N(ii)=length(t{ii});
        idx=find(diff(t{ii})==0);
        arcposition{ii}=[1 idx+1;idx N(ii)];
    end
end

arcpos=arcposition;


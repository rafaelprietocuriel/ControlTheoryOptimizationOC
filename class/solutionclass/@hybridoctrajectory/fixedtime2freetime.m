function hocTrj=fixedtime2freetime(hocTrj,jumppos,jumpid)

jumparg=jumpargument(hocTrj);
if nargin==1
    jumppos=find(jumparg<0);
    jumpid=[];
end
if nargin==2
    jumpid=[];
end
if isempty(jumppos)
    return
end

for ii=1:length(jumppos)
    if jumppos(ii)==length(jumparg)
        hocTrj.x=[hocTrj.x hocTrj.x([end end])+1];
        hocTrj.y=[hocTrj.y hocTrj.y(:,[end end])];
        hocTrj.arcinterval=[hocTrj.arcinterval hocTrj.arcinterval(end)];
        hocTrj.arcarg=[hocTrj.arcarg hocTrj.arcarg(end)];
        if isempty(jumpid)
            hocTrj.jumparg(jumppos(ii))=abs(hocTrj.jumparg(jumppos(ii)));
        else
            hocTrj.jumparg(jumppos(ii))=jumpid(ii);
        end
        hocTrj.jumparg=[hocTrj.jumparg 0];
        hocTrj.arcposition=[hocTrj.arcposition hocTrj.arcposition(2,end)+[1;2]];
    elseif jumppos(ii)==1
        hocTrj.x=[hocTrj.x([1 1]) hocTrj.x+1 ];
        hocTrj.y=[hocTrj.y(:,[1 1]) hocTrj.y];
        hocTrj.arcinterval=[hocTrj.arcinterval(1) hocTrj.arcinterval];
        hocTrj.arcarg=[ hocTrj.arcarg(1) hocTrj.arcarg];
        if isempty(jumpid)
            hocTrj.jumparg(jumppos(ii))=abs(hocTrj.jumparg(jumppos(ii)));
        else
            hocTrj.jumparg(jumppos(ii))=jumpid(ii);
        end
        hocTrj.jumparg=[0 hocTrj.jumparg];
        hocTrj.arcposition=[[1;2] hocTrj.arcposition+2];

    else
        if isempty(jumpid)
            hocTrj.jumparg(jumppos(ii))=abs(hocTrj.jumparg(jumppos(ii)));
        else
            hocTrj.jumparg(jumppos(ii))=jumpid(ii);
        end
    end
end
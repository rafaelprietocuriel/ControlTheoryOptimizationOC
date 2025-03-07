function initline=newsectionline(fid_read)
counter=0;
maxcounter=50;
if feof(fid_read)~=0
    initline='';
    return
end
while 1
    counter=counter+1;
    if counter>maxcounter
        ocmaterror('')
    end
    initline=strtrim(fgetl(fid_read));
    fidx=findstr(initline,'%');
    if ~isempty(fidx)
        initline(fidx(1):end)=[];
    end
    initline=strtrim(initline);
    if endsection(initline)
        initline='';
        break
    end
    if ~isempty(initline)
        break;
    end
    if feof(fid_read)
        initline='';
        return
    end
end

function b=endsection(inputstr)
fidx=findstr(inputstr,'%');
if ~isempty(fidx)
    inputstr(fidx(1):end)=[];
end
inputstr=strtrim(inputstr);
b=strcmp(inputstr,'\\');
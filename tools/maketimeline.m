function maketimeline(timelinefile,fignumber,stops,varargin)

filedirectory='';
framerate=[];
if nargin>=4
    filedirectory=varargin{1};
end
if nargin>=5
    framerate=varargin{12};
end

if isempty(filedirectory)
    filedirectory='.\';
end
if isempty(framerate)
    framerate=15;
end

fid=fopen(fullfile(filedirectory,[timelinefile '.txt']),'w');

counter=0;
stopflag=0;
for ii=1:fignumber
    if ii==1 || stopflag
        framerateentry=framerate;
    else
        framerateentry=[];
    end
    if any(ii==stops)
        fprintf(fid,'*:%d:%d\n',framerateentry,counter);
        stopflag=1;
    elseif stopflag
        fprintf(fid,':%d:%d\n',framerateentry,counter);
        stopflag=0;
    else
        fprintf(fid,':%d:%d\n',framerateentry,counter);
    end
    counter=counter+1;
end

fclose(fid);
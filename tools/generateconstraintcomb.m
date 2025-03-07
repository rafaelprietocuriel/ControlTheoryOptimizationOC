function str=generateconstraintcomb(combstr,outflag,varargin)
%
% STR=GENERATECONSTRAINTCOMB({'CC1','CC2','CC3'},1);
%
% STR=GENERATECONSTRAINTCOMB({'CC1','CC2','CC3'},1,'include',[1 2]); include
% constraints 1 and 2.
%
% STR=GENERATECONSTRAINTCOMB({'CC1','CC2','CC3','CC4'},1,'exclude',[1 3]); include
% constraints 1 and 2.
startidx=[];
include=[];
exclude=[];
if nargin==1
    outflag=1;
end
for ii=1:2:length(varargin)
    if ischar(varargin{ii})
        eval([varargin{ii} '=varargin{ii+1};'])
    end
end
if isempty(startidx)
    startidx=0;
end
comb=1:length(combstr);

totnum=length(combstr);

ctr=0;
for ii=1:totnum
    tmpcomb=nchoosek(comb,ii);
    for jj=1:size(tmpcomb,1)
        detect=0;
        if ~isempty(include)
            ctr0=0;
            while 1
                ctr0=ctr0+1;
                if any(strcmp({combstr{tmpcomb(jj,:)}},combstr(include(ctr0))))
                    detect=1;
                    break
                end
                if ctr0+1>length(include)
                    break
                end
            end

        else
            detect=1;
        end
        if ~isempty(exclude) && detect
            ctr0=0;
            while 1
                ctr0=ctr0+1;
                logicflag=any(strcmp({combstr{tmpcomb(jj,:)}},combstr(exclude(ctr0,1))));
                for kk=2:size(exclude,2)
                    logicflag=logicflag && any(strcmp({combstr{tmpcomb(jj,:)}},combstr(exclude(ctr0,kk))));
                end
                if logicflag
                    detect=0;
                    break
                end
                if ctr0+1>size(exclude,1)
                    break
                end
            end
        end
        if detect
            ctr=ctr+1;
            str{ctr}=[num2str(ctr+startidx) '::'];
            for kk=1:size(tmpcomb,2)
                if kk<size(tmpcomb,2)
                    str{ctr}=[str{ctr} combstr{tmpcomb(jj,kk)},','];
                else
                    str{ctr}=[str{ctr} combstr{tmpcomb(jj,kk)}];
                end
            end
        end
    end
end
if startidx==0
    str=[{'0::[]'} str];
end

if outflag
    for ii=1:length(str)
        fprintf('%s\n',str{ii})
    end
end
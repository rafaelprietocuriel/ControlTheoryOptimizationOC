function varargout=plotarrow(lhdl,pos,or,relarlegth)

%lhdl: line handle
%relarlegth: in percentage of curve length
%pos: in relative length of curve length: value betwenn 0 and 1 
%or: orientation positive or negative
if isempty(lhdl)
    if nargout==1
        varargout{1}=[];
    end
    return
end

if ~(ishandle(lhdl) && strcmp(get(lhdl,'Type'),'line'))
    if nargout==1
        varargout{1}=[];
    end
    return
end

if isempty(or)
    or=1;
end
ah=get(lhdl,'Parent');
ldata=[get(lhdl,'XData'); ...
    get(lhdl,'YData'); ...
    get(lhdl,'ZData')];
lim=[get(ah,'XLim'); ...
    get(ah,'YLim'); ...
    get(ah,'ZLim')];
Dlim=diff(lim,[],2);

ldata(:,ldata(1,:)<lim(1,1))=[];
ldata(:,ldata(2,:)<lim(2,1))=[];
ldata(:,ldata(1,:)>lim(1,2))=[];
ldata(:,ldata(2,:)>lim(2,2))=[];
if size(ldata,1)==3
    fprintf('Function properly only works for 2D plots\n')
    if nargout==1
        varargout{1}=[];
    end
    return
end

oldunits = get(ah, 'Units');
set(ah, 'Units', 'Normalized');
axpos=get(ah,'Position');
set(ah, 'Units', oldunits);
ax_per_data=axpos(3:4).';
ax_per_data=ax_per_data./Dlim(1:2);
axpos=axpos(1:2).';

if isempty(ldata)
    if nargout==1
        varargout{1}=[];
    end
    return
end
lgth=[0 cumsum(sqrt(sum(diff(ldata,[],2).^2)))];
[lgth,idx]=unique(lgth);
ldata=ldata(:,idx);
abspos=pos*lgth(end);
absarlegth=relarlegth*lgth(end)/100;

data=interp1(lgth,ldata.',[abspos-absarlegth/2 abspos+absarlegth/2]).';
ndata=[axpos axpos]+(data-[lim(1:2,1) lim(1:2,1)]).*([ax_per_data ax_per_data]);

if or<0
    ndata=ndata(:,[2 1]);
end

ar=annotation('arrow',ndata(1,:),ndata(2,:));
set(ar,'Color',get(lhdl,'Color'))

if nargout==1
    varargout{1}=ar;
end

%for testing
%  hold on
%  plot(ah,ldata(1,[1 end]),ldata(2,[1 end]),'xr');
%  plot(ah,data(1,:),data(2,:),'xk');

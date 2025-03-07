function savefigure(fighdl,figfilename,imgtype,varargin)
%
% SAVEFIGURE prints a figure to an image.
%
% SAVEFIGURE(H) prints a figure plotted in figure, with handle H. It is the
% current figure by default.
%
% SAVEFIGURE(H,FN)
%
% SAVEFIGURE(H,FN,IMGT)
%

if nargin==2
    imgtype='pdf';
end
if isempty(figfilename)
    ocmaterror('No file name provided.')
end
if isempty(fighdl)
    fighdl=gcf;
end
if isempty(imgtype)
    imgtype='pdf';
end
if ~ishandle(fighdl)
    ocmaterror('First argument is not a valid graphic handle.')
end
if ~strcmpi(get(fighdl,'Type'),'figure')
    ocmaterror('First argument is not a valid figure handle.')
end
if ~any(strncmp(varargin,'-r',2))
    varargin{end+1}='-r750';
end
varargin(strncmp(varargin,'-d',2))=[];

% if 'png' set '-zbuffer '
if ~verLessThan('Matlab','8.4')
    fighdl=fighdl.Number;
end

if strcmp(imgtype,'fig')
    saveas(fighdl,figfilename)
else
    print(['-d' imgtype],['-f' num2str(fighdl)],varargin{:},figfilename)
end

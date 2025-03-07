function savefigure(odeObj,varargin)
%
% SAVEFIGURE prints a figure to an image.
%
% SAVEFIGURE(H) prints a figure plotted in figure, with handle H. It is the
% current figure by default.
%
% SAVEFIGURE(H,IMGT) possible image types are 'bmp', 'eps' and 'png'. The
% default value is 'eps'.
%
% SAVEFIGURE(H,IMGT,FN,FF) with FN a file name and with FF a folder can be
% specified. The defualt folder is the data folder of OCOBJ.
%
% SAVEFIGURE(H,IMGT,FN,FF,RES) with RES the resolution of the image can be
% specified. By the default RES = 750.


figfilename='';
imgtype='';
fighdl=[];
resolution=[];
savedir='';
if nargin>=2
   fighdl=varargin{1};
end
if nargin>=3
   imgtype=varargin{2};
end
if nargin>=4
   figfilename=varargin{3};
end
if nargin>=5
   savedir=varargin{4};
end
if nargin>=6
   resolution=varargin{5};
end
if isempty(savedir)
    savedir=fullocmatfile(userdatafolder(odeObj)); % Data directory
end
if isempty(figfilename)
    figfilename=filename(odeObj);
end
if isempty(figfilename)
    ocmaterror('For empty oc model filename for figure has to be provided.');
end
if isempty(fighdl)
    fighdl=gcf;
end
if isempty(imgtype)
    imgtype='pdf';
end
if isempty(resolution)
    resolution=750;
end
% if 'png' set '-zbuffer '

savefigure(fighdl,fullfile(savedir,[figfilename '.' imgtype]),imgtype,['-r' num2str(resolution)],varargin{6:end})
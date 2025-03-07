function out = setocmatoptions(varargin)

if mod(nargin,2) 
    error('Number of input arguments must be in form of (''Fieldname'',Value) or empty.')
end


load OCMatOptions
mainoptionsdir=getocmatpath('options');

if nargin == 0
    out=ocmatoptions;
    return
end    

for ii = 1:2:nargin
    fieldname=varargin{ii};
    value=varargin{ii+1};
    if fieldcontrol(fieldname,value)
        ocmatoptions=setfield(ocmatoptions,fieldname,value);
    end
end

out=ocmatoptions;
save(fullfile(mainoptionsdir,['OCMatOptions.mat']),'ocmatoptions')


   


function out = getocmatoptions(varargin)

idx=1;

load OCMatOptions


if nargin == 0
    out=ocmatoptions;
    return
end    

for ii = 1:nargin
    fieldname=varargin{ii};
    
    if fieldcontrol(fieldname)
        out{idx}=getfield(ocmatoptions,fieldname);
        idx=idx+1;
    end
end



   


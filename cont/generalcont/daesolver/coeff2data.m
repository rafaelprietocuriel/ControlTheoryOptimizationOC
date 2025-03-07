function [x,y,z1,freepar]=coeff2data(tmesh,coeff)
global OCMATCONT
persistent continuationindex freeparameterindex ycoefficientidx firstorderzcoefficientidx

if continuationindex || continuationindex~=length(coeff)
    ycoefficientidx
    firstorderzcoefficientidx
end
y=coeff(OCMATCONT.MESHDATA.ycoefficientidx); % y coefficients
z1=coeff(OCMATCONT.MESHDATA.firstorderzcoefficientidx); % z coefficients for the first order components

x=coeff2points_new(tmesh,coeff,'collocationgrid');
freepar=coeff(end-OCMATCONT.freeparameternum+1:end);
function out=freeparameters(ppdeTrj)
out=[];

if isfield(ppdeTrj.discretizationinfo,'freeparameters')
    out=ppdeTrj.discretizationinfo.freeparameters(:).';
end
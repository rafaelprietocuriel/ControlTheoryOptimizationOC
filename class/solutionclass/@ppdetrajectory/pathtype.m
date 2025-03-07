function out=pathtype(ppdeTrj)

out=[];
if isfield(ppdeTrj.discretizationinfo,'pathtype')
    out=ppdeTrj.discretizationinfo.pathtype;
end
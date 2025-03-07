function out=initialtime(ocgTrj)

if ~isfield(ocgTrj.initargument,'t')
    out=[];
else
    out=ocgTrj.initargument.t;
end
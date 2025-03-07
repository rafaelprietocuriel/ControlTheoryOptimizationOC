function out=independentvar(docAsym)

out=independentvar(docAsym.doctrajectory);
if strcmp(pathtype(docAsym),'u')
    out=out(end:-1:1);
end


function out=solverinfo(ocgTrj)

deg=multiplicity(ocgTrj);
if deg>1
    out=cell(1,deg);
    for ii=1:deg
        out{ii}=ocgTrj(ii).solverinfo;
    end
else
    out=ocgTrj.solverinfo;
end
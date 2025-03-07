function b=isocsolutionclass(a)

b=1;
folder=dir('./class/solutionclass/@*');

for ii=1:numel(folder)
    b=b||isa(a,folder(ii).name(2:end));
end
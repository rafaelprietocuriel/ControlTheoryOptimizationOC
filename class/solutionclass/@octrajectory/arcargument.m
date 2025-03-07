function out=arcargument(ocTrj)

if strcmp(solver(ocTrj),'dae')
    out=1;
else
    out=ocTrj.arcarg;
end
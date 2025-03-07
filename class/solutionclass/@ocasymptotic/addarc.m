function ocAsym=addarc(ocAsym,arc,iposition)

if nargin==2
    iposition=1;
end

ocAsym=ocasymptotic(addarc(octrajectory(ocAsym),arc,iposition),ocAsym.dynprimitive);
function ocTrj=addarc(ocTrj,arc,iposition)

if nargin==2
    iposition=1;
end

ocTrj=octrajectory(addarc(struct(ocTrj),arc,iposition));
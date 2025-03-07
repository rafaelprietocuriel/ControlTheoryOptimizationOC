function out=arcinterval(ocTrj,idx)

if nargin==1
    idx=1:length(ocTrj.arcinterval);
end
out=ocTrj.arcinterval(idx);
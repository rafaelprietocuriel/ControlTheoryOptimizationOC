function out=fundamentalsolution(ocTrj)
out=[];
if isempty(findstr(conttype(ocTrj),'F'))
    return
end
sinfo=solverinfo(ocTrj);
out=ocTrj.y(sinfo.fundcoord,:);
out=reshape(out,sqrt(sinfo.fundnumcoord),sqrt(sinfo.fundnumcoord),[]);
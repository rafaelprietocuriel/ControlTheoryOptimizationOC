function ocTrj=addobjectivevaluepath(ocTrj,objpath,objcoord)

if isempty(ocTrj)
    return
end

ocTrj.y(objcoord,:)=objpath;
ocTrj.solverinfo.objecivevaluecoordinate=objcoord;

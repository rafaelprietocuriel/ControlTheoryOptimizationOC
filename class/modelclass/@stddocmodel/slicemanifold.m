function out=slicemanifold(docObj,idx)

ocRes=contresult(docObj);
if isempty(idx)
    out=[];
    return
end
if ~isnumeric(idx) || numel(idx)~=1 || idx<1
    ocmaterror('')
end
if idx>length(ocRes)
    ocmaterror('')
end

switch ocRes{idx}.ContinuationClassification
    case {'extremal2fp','dextremal2fp'}
        coord=[statecoord(docObj).';costatecoord(docObj).'];
        numcoord=length(coord);
        numidx=1:numcoord;
        numcontsol=length(ocRes{idx}.ContinuationSolution);
        smf.y=zeros(numcoord+1,numcontsol);
        smf.userinfo.continuationparameter=zeros(1,numcontsol);
        smf.arcarg=zeros(1,numcontsol);
        for ii=1:numcontsol
            smf.y(numidx,ii)=ocRes{idx}.ContinuationSolution(ii).y0(coord,1);
            smf.y(numcoord+1,ii)=calcobjectivevalue(docObj,doctrajectory(ocRes{idx}.ContinuationSolution(ii)));
            smf.userinfo.continuationparameter(ii)=ocRes{idx}.ContinuationSolution(ii).solverinfo.continuationparameter;
            smf.arcarg(ii)=ocRes{idx}.ContinuationSolution(ii).arcarg(1);
        end
        smf.arcposition=[1;size(smf.y,1)];
        out=occurve(smf);
end
function ocLC=homoclinic2lc(ocMP,varargin)

ocLC=[];
if ~isocmultipath(ocMP)
    return
end
if multiplicity(ocMP)>2
    return
end

% find unstable path
uidx=[];
ocEx=cell(2,1);
for ii=1:multiplicity(ocMP)
    ocEx{ii}=octrajectory(ocMP.solutionclass{ii});
    if strcmp(pathtype(ocMP.solutionclass{ii}),'u')
        uidx=ii;
    end
end

if isempty(uidx)
    ocmatmsg('ocmultipath is not a homoclinic path.')
    return
end
    
sidx=1:2;
sidx(uidx)=[];

ocEx=ocEx([sidx uidx]);
ocEx{2}=reversetime(ocEx{2});

ocLC.octrajectory=horzcat(ocEx{:});
ocLC.octrajectory=mergearc(ocLC.octrajectory,arcnum(ocEx{1}));
ocLC.period=ocLC.octrajectory.arcinterval(end);
ocLC=dynprimitive(ocLC);
function ocTrj=ocmultipath2octrajectory(ocMP,varargin)

ocTrj=octrajectory([]);
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
    return
end
    
sidx=1:2;
sidx(uidx)=[];

ocEx=ocEx([sidx uidx]);
ocEx{2}=reversetime(ocEx{2});

ocTrj=horzcat(ocEx{:});
ocTrj=mergearc(ocTrj,arcnum(ocEx{1}));
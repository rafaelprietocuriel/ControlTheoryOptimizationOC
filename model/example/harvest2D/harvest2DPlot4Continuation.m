function out = harvest2DPlot4Continuation(t,dynVar,tangent,contpar,flag)
%
% the BC file returning the continuation residuum
                                                                                                                                                                                                      
% global variable
global OcBVPVar OcSolverVar AsymVar LAsymVar %LAsymVar2
out=[];
switch OcBVPVar.problem
    case 'asymptoticextremal'
        try
            AsymVar.initpt=[AsymVar.initpt dynVar(1:4,1)];
        catch
            AsymVar.initpt=dynVar(1:4,1);
        end
        switch numel(OcBVPVar.ArcIdentifier)
            case 1
                subplot(2,2,1)
                plot(AsymVar.initpt(1,:),AsymVar.initpt(3,:)),
                subplot(2,2,2)
                plot(AsymVar.initpt(2,:),AsymVar.initpt(3,:)),
                subplot(2,2,3)
                plot(AsymVar.initpt(1,:),AsymVar.initpt(2,:)),
                subplot(2,2,4)
                plot(dynVar(1,:),dynVar(2,:)),
        end
    case 'limitpointextremal'
        LAsymVar.initpt=[LAsymVar.initpt dynVar(1:4,1)];
        switch numel(OcBVPVar.ArcIdentifier)
            case 1
                dynVar=reshape(dynVar,4,OcSolverVar.meshnum);
                subplot(1,2,1)
                plot(LAsymVar.initpt(1,:),LAsymVar.initpt(2,:))
                subplot(1,2,2)
                plot(dynVar(1,:),dynVar(2,:)),
        end
end
drawnow
figure(gcf)

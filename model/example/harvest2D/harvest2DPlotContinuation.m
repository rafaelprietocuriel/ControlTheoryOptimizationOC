function h = harvest2DPlotContinuation(t,dynVar,pararg,arcarg,freepar,tangent,plotflag)
global OCMATCONT
h=[];
leftarcindex=OCMATCONT.HE.TIMEDDATA.leftarcindex;
rightarcindex=OCMATCONT.HE.TIMEDDATA.rightarcindex;

numarc=numel(arcarg);
if numarc==1
    h=plot(dynVar(1,:),dynVar(2,:));
elseif numarc==2
    idx1=leftarcindex(1):rightarcindex(1);
    idx2=leftarcindex(2):rightarcindex(2);
    if arcarg(1)==1
        h(1)=plot(dynVar(1,idx1),dynVar(2,idx1),'b');
        hold on
        h(2)=plot(dynVar(1,idx2),dynVar(2,idx2),'r');
        hold off
    else
        h(1)=plot(dynVar(1,idx1),dynVar(2,idx1),'r');
        hold on
        h(2)=plot(dynVar(1,idx2),dynVar(2,idx2),'b');
        hold off
    end
end
xlabel('$F$','Interpreter','Latex')
ylabel('$A$','Interpreter','Latex')
drawnow
figure(gcf)

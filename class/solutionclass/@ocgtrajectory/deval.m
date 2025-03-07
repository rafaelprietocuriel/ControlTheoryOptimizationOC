function ocTrj=deval(ocTrj,xint)
%
% DEVAL evaluate an octrajectory
%
% DEVAL(OCTRJ,XINT) evaluates the instance OCTRJ of an octrajectory at the
% grid XINT. XINT has to consist at least of the initial and endtime of
% OCTRJ. If OCTRJ consists of multiple arcs XINIT has to contain also the
% grid points of the switching points between the different arcs.
%
% OCTRJN=DEVAL(OCTRJ,XINT) OCTRJN is an instance of octrajectory evaluated
% at XINT. The 'solver' field is set to empty.

sol=octrajectory2nativematlab(ocTrj);
sol.solver='gbvp4c';

ocTrjStruct=struct(ocTrj);
odenum=ocTrjStruct.odenum;
ocTrjStruct=struct(ocTrjStruct.octrajectory);
ocTrjStruct.odenum=odenum;
ocTrjStruct.solver='';
arcpos=arcposition(ocTrj);
% if size(arcpos,2)>1
%     for ii=2:size(arcpos,2)
%         xint=[xint ocTrjStruct.x([arcpos(2,ii-1) arcpos(1,ii)])];
%     end
%     xint=sort(xint);
% end
arcpos2=find(diff(xint)==0);
arcpos2=[[1 arcpos2+1];[arcpos2 numel(xint)]];
timehorizon=ocTrjStruct.timehorizon;
switch sol.solver
    case 'gbvp4c'
        for ii=1:size(arcpos,2)
            Sxint=interp1(sol.x(arcpos(1,ii):arcpos(2,ii)),sol.y{ii}.',xint(arcpos2(1,ii):arcpos2(2,ii))).';
            ocTrjStruct.y{ii}=Sxint;
        end
end
ocTrjStruct.arcposition=arcpos2;
ocTrjStruct.timehorizon=timehorizon;
ocTrjStruct.x=xint;
ocTrjStruct.arcinterval(end)=timehorizon;
ocTrj=ocgtrajectory(ocTrjStruct);

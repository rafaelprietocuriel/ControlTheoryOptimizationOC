function ocAsym=deval(ocAsym,xint)

ocAsym.ocgtrajectory=deval(ocAsym.ocgtrajectory,xint(xint<=1));
ocAsym.ocgtrajectory.arcinterval(end)=ocAsym.ocgtrajectory.timehorizon*xint(end);
ocAsym.ocgtrajectory.timehorizon=ocAsym.ocgtrajectory.timehorizon*xint(end);
ocAsym.ocgtrajectory.x=xint/xint(end);
xint(xint<=1)=[];
ocAsym.ocgtrajectory.y=[ocAsym.ocgtrajectory.y repmat(ocAsym.limitset.y,1,length(xint))];
ocAsym.ocgtrajectory.arcposition(2,end)=ocAsym.ocgtrajectory.arcposition(2,end)+length(xint);

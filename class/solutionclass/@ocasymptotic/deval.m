function ocAsym=deval(ocAsym,xint)

if max(xint)<=ocAsym.octrajectory.x(end)
    ocAsym.octrajectory=deval(ocAsym.octrajectory,xint);
else
    if size(ocAsym.limitset.y,2)<2
        ocAsym.octrajectory=deval(ocAsym.octrajectory,xint(xint<=1));
        ocAsym.octrajectory.arcinterval(end)=ocAsym.octrajectory.timehorizon*xint(end);
        ocAsym.octrajectory.timehorizon=ocAsym.octrajectory.timehorizon*xint(end);
        ocAsym.octrajectory.x=xint/xint(end);
        xint(xint<=1)=[];
        ocAsym.octrajectory.y=[ocAsym.octrajectory.y repmat(ocAsym.limitset.y,1,length(xint))];
        ocAsym.octrajectory.arcposition(2,end)=ocAsym.octrajectory.arcposition(2,end)+length(xint);
    else
    end
end

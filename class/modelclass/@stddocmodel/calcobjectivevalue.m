function out=calcobjectivevalue(docObj,docTrj)

objval=objectivevalue(docObj,docTrj,1);
out=sum(objval);
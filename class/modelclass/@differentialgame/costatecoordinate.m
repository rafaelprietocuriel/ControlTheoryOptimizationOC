function out=costatecoordinate(dgObj,player)
% player is the number of a "real" player of the differential game, i.e.
% player 0 is not allowed
if nargin==1
    player=[];
end

n=playernum(dgObj);
statenum=retrievedifferentialgameinformation(dgObj.Model,'statenum');
statenum=sum(statenum.value);
costatenum=retrievedifferentialgameinformation(dgObj.Model,'costatenum');
costatenum=sum(costatenum.value);
if isempty(dgObj)
    out=[];
    return
end

costatecoordinate=statenum+reshape(1:costatenum,[],n);
if isempty(player)
    out=costatecoordinate(1,1):costatecoordinate(end,n);
else
    out=costatecoordinate(:,player);
end
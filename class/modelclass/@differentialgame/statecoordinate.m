function out=statecoordinate(dgObj,player)

if nargin==1
    player=[];
end
out=[];
if isempty(dgObj)
    return
end
existplayer0=retrievedifferentialgameinformation(dgObj.Model,'existplayer0');
statenum=retrievedifferentialgameinformation(dgObj.Model,'statenum');
statenum=[0 cumsum(statenum.value)];
if existplayer0.value
    player=player+1;
end
if isempty(player)
    out=1:sum(statenum(end));
else
    out=statenum(player)+1:statenum(player+1);
end
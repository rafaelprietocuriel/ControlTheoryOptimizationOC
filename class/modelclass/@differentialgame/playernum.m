function num=playernum(dgObj)
%
% number of players

num=[];
if isempty(dgObj)
    return
end
info=retrievedifferentialgameinformation(dgObj.Model,'playernum');
num=info.value;
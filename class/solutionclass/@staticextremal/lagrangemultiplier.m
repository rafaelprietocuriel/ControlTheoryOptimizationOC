function out=lagrangemultiplier(statEx)
out=[];
if isempty(statEx)
    return
end
out=statEx.lm;

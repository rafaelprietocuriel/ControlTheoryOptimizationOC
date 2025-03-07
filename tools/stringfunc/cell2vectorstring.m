function string=cell2vectorstring(cellstring)

if isempty(cellstring)
    string=cellstring;
    return
elseif ischar(cellstring)
    string=['[' cellstring ']'];
    return
end
string=['[' cellstring{1}];
for ii=2:numel(cellstring)
    if ~isempty(cellstring{ii})
        string=[string ',' cellstring{ii}];
    end
end
string=[string ']'];

function string=removesquarebracket(string)
% the square brackets are removed
string=regexprep(string,'\[||\]','');
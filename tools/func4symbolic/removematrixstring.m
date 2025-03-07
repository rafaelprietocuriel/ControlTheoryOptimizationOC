function string=removematrixstring(string)
% the signal word 'matrix' together with the outer brackets [] are removed
% from strings returned by the symbolic toolbox
string=regexprep(string,'\<(matrix\(\[)||(\]\))\>','');
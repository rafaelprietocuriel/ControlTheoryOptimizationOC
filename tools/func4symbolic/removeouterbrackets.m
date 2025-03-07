function s = removeouterbrackets(s)
%TRIM  TRIM(s) deletes any white spaces.
s=regexprep(s,'(^[)|(]$)','');


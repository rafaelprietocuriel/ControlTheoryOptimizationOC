function [matchname,leftidx,rightidx]=findvariableintext(filetext,variabletype)

if nargin==1
    variabletype='OCMATVARIABLE';
end

switch variabletype
    case 'OCMATVARIABLE'
        [matchname,leftidx,rightidx]=regexp(filetext,'\$[A-Z0-9]+\$','match');
        
    case 'FORNEXTCOUNTER'
        [matchname,leftidx,rightidx]=regexp(strrep(filetext,'!FOR',''),'=','split');
        matchname(2)=[];
        matchname=strtrim(matchname);
        
    case 'EVALCOMMAND'
        [matchname,leftidx,rightidx]=regexp(filetext,'\!EVAL(.+)\!','match');

    case 'IFCLAUSE'
        [matchname,leftidx,rightidx]=regexp(filetext,'\!(ENDIF)*(IF\ [\~0-9<>-+/=*()&\ ]+)*\!','match');
        
    case 'ADDARCCASE'
        [matchname,leftidx,rightidx]=regexp(filetext,'\!(START)*(END)*ADDARCCASE\!','match');
        
    case 'ADDCASE'
        [matchname,leftidx,rightidx]=regexp(filetext,'\!(START)*(END)*ADDCASE\!','match');
        
    case 'FORNEXT'
        [matchname,leftidx,rightidx]=regexp(filetext,'\!(ENDFOR)*(FOR\ [\w =:\[\]]*)*\!','match');
end
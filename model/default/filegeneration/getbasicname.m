function basename=getbasicname(basename)
%
% returns the variable names which are used for the auomatically generated
% files
switch basename
    case 'independent'
        basename='t';
        
    case 'time'
        basename='t';
        
    case 'connectiontime'
        basename='ct';
        
    case 'space'
        basename='z';
        
    case 'spacemid'
        defaultspacearg=getbasicname('space');
        basename=[defaultspacearg 'mid'];
        
    case 'state'
        basename='x';
        
    case 'costate'
        basename='lambda';
        
    case 'dependent'
        % in general the dependent variable consists of state, costate and
        % implicit control/Lagrangemultipliers
        basename='depvar';
        
    case 'control'
        basename='ctrl';
        
    case 'controldx'
        basename='dctrldx';

    case 'lagrangemultcc'
        basename='lagmcc';

    case 'constraint'
        basename='cstr';

    case 'lagrangemultsc'
        basename='lagmsc';

    case 'parametervariables'
        basename='par';

    case 'arcidentifiervar'
        basename='arcid';

    case 'jumpidentifiervar'
        basename='jumpid';
        
    case 'endtime'
        basename='T';
        
    case 'impulsetime'
        basename='ttau';
        
    case 'discountrate'
        basename='r';
        
    case 'exogenousfunction'
        basename='exfunc';
        
    case 'exogenousstate'
        basename='depvar';
        
    case 'nonsmoothfunction'
        basename='nonsfunc';
        
    case 'femdata'
        basename='femdata';
        
    case 'femdatagrid'
        femdata=getbasicname('femdata');
        basename=[femdata '.grid'];
        
    case 'femdatagridmid'
        femdata=getbasicname('femdata');
        basename=[femdata '.gridmid'];

    case 'femdatadgrid'
        femdata=getbasicname('femdata');
        basename=[femdata '.dgrid'];

    case 'femdatagridnum'
        femdata=getbasicname('femdata');
        basename=[femdata '.gridnum'];

    case 'femdatagridnumm1'
        femdata=getbasicname('femdata');
        basename=[femdata '.gridnumm1'];

    case 'femdataxcoordm1'
        femdata=getbasicname('femdata');
        basename=[femdata '.xcoordm1'];

    case 'femdataxcoordp1'
        femdata=getbasicname('femdata');
        basename=[femdata '.xcoordp1'];

    case 'femdatainvMK'
        femdata=getbasicname('femdata');
        basename=[femdata '.invMK'];

    case 'impulsedependent'
        % the impulsedependent variable consists of the left and right side
        % limits at the impulse time for state and costates 
        basename='idepvar';
        
    case 'impulsecontrol'
        basename='ictrl';
        
    case 'locstate'
        basename='loc_state';
        
    case 'loccostate'
        basename='loc_costate';
        
    case 'loccontrol'
        basename='loc_control';
        
    case 'gradientloccontrol'
        basename='gradloc_control';

    case 'lagrangemult'
        basename='lm';

    case 'lagrangemultineq'
        basename='lm';

    case 'lagrangemulteq'
        basename='lmeq';
                
    case 'variable'
        basename='x';

end
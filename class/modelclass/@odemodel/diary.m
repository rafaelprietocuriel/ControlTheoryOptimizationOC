function diary(odeObj,flag)
%
% DIARY toggles model diary file on/off
%
% DIARY(OCOBJ) writes to the model specific diary file
% '[modelname]Diary.m' in the model data folder. New commands are appended
% to an already existing diary file. Calling DIARY(OCOBJ) a second time
% turns the diary modus off.

if isempty(odeObj)
    ocmatmsg('Oc model is empty.\n')
    return
end

if nargin==1
    switch get(0,'Diary')
        case 'off'
            flag='on';
        case 'on'
            flag='off';
    end
end

modeldiaryfile=fullocmatfile(userdatafolder(odeObj),[modelname(odeObj) 'Diary.m']);
switch flag
    case 'on'
        diary(modeldiaryfile)
    case 'off'
        actdiaryfile=get(0,'DiaryFile');
        if ~strcmp(modeldiaryfile,actdiaryfile)
            ocmatmsg('Actual diary file ''%s'' is not the standard model diary file\n',actdiaryfile)
            answer=input('Interrupt?  y/(n): ','s');
            if isempty(answer)
                % default value 'y'
                answer='n';
            end
            if strcmpi(answer,'y')
                disp('Diary not set to off.')
                return
            end
        end
        diary off
    otherwise
        ocmatmsg('Input argument ''%s'' is unknown.',flag)
end
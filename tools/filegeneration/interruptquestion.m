function b=interruptquestion()
%
%
while 1
    accept=lower(input('Interrupt: y/(n) : ','s'));
    if isempty(accept)
        accept='n';
    end
    if strcmpi(accept,'y')
        b=1;
        break
    elseif strcmp(accept,'n')
        b=0;
        break
    end
end

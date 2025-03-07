function [symkernel,symbolicinfo]=getsymkernel()

% determine if symbolic toolbox is registered
symbolicinfo=ver('symbolic');
if ~isempty(symbolicinfo)
    if verLessThan('matlab','7.6.0')
        % maple kernel is used for symbolic calculations
        symkernel='maple';
    else
        % mupad kernel is used for symbolic calculations
        symkernel='mupad';
    end
else
    symkernel='';
end
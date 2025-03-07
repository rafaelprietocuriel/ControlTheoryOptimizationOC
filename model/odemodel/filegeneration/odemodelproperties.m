function [mainsectionnames mandatorysectionnames]=odemodelproperties(mainsection)
%
%
mandatorysectionnames='';
if nargin==0
    mainsection='init';
end

switch mainsection
    case 'init'
        mainsectionnames={'description','variable','dynamics','parameter','exogenousfunction'};
        mandatorysectionnames=mainsectionnames(1:3);
end

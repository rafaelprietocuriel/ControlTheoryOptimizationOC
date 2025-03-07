function extremal=generategradextremal(ocgTrj,varargin)
% generate extremal for gradient initialization

extremal=[];
if nargin==1
    normalizetime=0;
else
    normalizetime=varargin{1};
end
if isempty(normalizetime)
    normalizetime=0;
end
if isempty(ocgTrj)
    return
end
var=variable(ocgTrj);

switch octype(ocgTrj)
    case 'static'
    otherwise
        arg=argument(ocgTrj);
        arg0=initargument(ocgTrj);

        fn=fieldnames(arg);
        for ii=1:length(fn)
            extremal.(fn{ii})=arg.(fn{ii});
        end
        fn=fieldnames(arg0);
        for ii=1:length(fn)
            extremal.(fn{ii})=arg0.(fn{ii});
        end
        extremal.timehorizon=timehorizon(ocgTrj);
        extremal.dv=[];
end
fn=fieldnames(var);
for ii=1:length(fn)
    extremal.(fn{ii})=var.(fn{ii});
end
extremal.objectivevalue=objectivevalue(ocgTrj);
if normalizetime
    if extremal.t(end)~=0
        extremal.t=extremal.t/extremal.t(end);
    end
end

switch octype(ocgTrj)
    case 'concentrated'
        extremal.coeff=[extremal.y(:);extremal.y(:,1);extremal.cst_y(:,end)];
    case 'static'
        %extremal.coeff=extremal.y(:);
end

function ocgTrj=deval(ocgTrj,tnew)

t=time(ocgTrj);

if numel(tnew)==1
    tnew=linspace(t(1),t(end),tnew);
end

if size(ocgTrj.variable.v,1)==1
    ocgTrj.variable.v=interp1(t,ocgTrj.variable.v,tnew);
else
    ocgTrj.variable.v=interp1(t,ocgTrj.variable.v.',tnew).';
end
if size(ocgTrj.variable.y,1)==1
    ocgTrj.variable.y=interp1(t,ocgTrj.variable.y,tnew);
else
    ocgTrj.variable.y=interp1(t,ocgTrj.variable.y.',tnew).';
end
if size(ocgTrj.variable.cst_y,1)==1
    ocgTrj.variable.cst_y=interp1(t,ocgTrj.variable.cst_y,tnew);
else
    ocgTrj.variable.cst_y=interp1(t,ocgTrj.variable.cst_y.',tnew).';
end

if isfield(ocgTrj.variable,'o')
    if size(ocgTrj.variable.o,1)==1
        ocgTrj.variable.o=interp1(t,ocgTrj.variable.o,tnew);
    else
        ocgTrj.variable.o=interp1(t,ocgTrj.variable.o.',tnew).';
    end
end
ocgTrj.argument.t=tnew;

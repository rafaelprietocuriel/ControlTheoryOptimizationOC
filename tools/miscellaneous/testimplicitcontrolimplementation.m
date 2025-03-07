function [dUdtsymori,dUdtsym,dUdtsym0,dUdX,dUdX0,JUDX,JUDX0]=testimplicitcontrolimplementation(arcid)

m=stdocmodel('nonlincapitalaccumulationi');
[dUdX,dUdtsym,JUDX]=generateimplicitcontroldynamics(m,arcid);

m0=stdocmodel('nonlincapitalaccumulation');
[dUdX0,dUdtsym0,JUDX0]=generateexplicitcontroldynamics(m0,arcid);

ctrlval=control(m0,[],arcid,1);
ctrl=retrievemodelinformation(m0,'controlname');
ctrlname=ctrl.value;
fs=ocmatfindsym(string2cell(removematrixstring(char(dUdtsym)),'matrix'),getsymkernel());
repflag=0;
dUdtsymori=dUdtsym;
for ii=1:length(ctrlval)
    if any(strcmp(fs,ctrlname{ii}))
        dUdtsym=subs(dUdtsym,sym(ctrlname{ii}),ctrlval(ii));
        repflag=1;
    end
end
if repflag
    fprintf('\nControl value was replaced in implicitly given control dynamics!\n')  
end
fprintf('\n\nTest for arcid = %d\n',arcid)

fprintf('\nControl dynamics:\n')
e=compare(dUdtsym,dUdtsym0);
switch e
    case 0
        fprintf('\tControl dynamics are equal.\n')
    case 1
        fprintf('\tControl dynamics differ in size.\n')
    case 2
        fprintf('\tControl dynamics have different entries.\n')
        tst=(dUdtsym-dUdtsym0);
        disp(tst)
        fprintf('\n')
    case 3
end

e=compare(dUdX,dUdX0);
fprintf('\nDUDX:\n')
switch e
    case 0
        fprintf('\tDUDX are equal.\n')
    case 1
        fprintf('\tDUDX differ in size.\n')
    case 2
        fprintf('\tDUDX have different entries.\n')
        tst=simple(dUdX-dUdX0);
        disp(tst)
        fprintf('\n')
    case 3
end

e=compare(JUDX,JUDX0);
fprintf('\nJacobian:\n')
switch e
    case 0
        fprintf('\tJacobian are equal.\n')
    case 1
        fprintf('\tJacobian differ in size.\n')
    case 2
        fprintf('\tJacobian have different entries.\n')
        tst=simple(JUDX-JUDX0);
        disp(tst)
        fprintf('\n')
    case 3
end




function e=compare(T1,T2)
e=0;
if any(size(T1)~=size(T2))
    e=1;
    return
end

DT=simple(T1-T2);
if ~isempty(findsym(DT))
    e=2;
    return
end

if max(double(DT(:)))
    e=3;
    return
end
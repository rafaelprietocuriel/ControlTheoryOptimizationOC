function ocAsym=periodicsolution2trajectory(ocLC,n)
%
%
arcarg=arcargument(ocLC);
arcint=arcinterval(ocLC);

numarc=length(arcarg);

t=independentvar(ocLC);
x=dependentvar(ocLC);
add2t=repmat((0:n-1).'*numarc,1,length(t)).';
add2t=add2t(:).';
ttot=repmat(t,1,n)+add2t;
xtot=repmat(x,1,n);
arcargtot=repmat(arcarg,1,n);
add2arcint=repmat((0:n-1).'*arcint(end),1,numarc).';
add2arcint=add2arcint(:).';
arcinttot=[0 repmat(arcint(2:end),1,n)+add2arcint];

arcposition=find(diff(ttot)==0);

ocTrj.x=ttot;
ocTrj.y=xtot;
ocTrj.arcarg=arcargtot;
ocTrj.arcinterval=arcinttot;
ocTrj.arcposition=[1 arcposition+1;arcposition length(ttot)];

mergidx=find(diff(arcargtot)==0);
ocTrj=octrajectory(ocTrj);

for ii=length(mergidx):-1:1
    ocTrj=mergearc(ocTrj,mergidx(ii));
end

ocAsym=ocasymptotic(ocTrj,ocLC);

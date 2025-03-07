function [b,val]=testjump(ocObj,hocTrj,varargin)
%
% TESTJUMP returns 1 if switching conditions at two adjacent arcs or jump
% conditions at the initial and end time are satisfied.

b=[];
val=[];
opt=[];
if isempty(ocObj)
   return
end
if ~ishybridoctrajectory(hocTrj)
    return
end
if nargin>2
    opt=varargin{1};
end
if isempty(opt)
    opt=defaultocoptions;
end
AdmissibleTolerance=getocoptions(opt,'GENERAL','AdmissibleTolerance');

par=parametervalue(ocObj);
coord=[statecoord(ocObj).';costatecoord(ocObj).'];
numj=jumpnum(hocTrj);
jumparg=jumpargument(hocTrj);
arcarg=arcargument(hocTrj);
b=ones(1,numj);
jdepvar=jumpdependentvar(hocTrj);
jt=arcinterval(hocTrj);
val=zeros(2*statenum(ocObj)+1,numj);
for ii=1:numj
    if jumparg(ii)
        if ii==1
            out1=feval(ocObj,'InteriorImpulseBC',jt(ii),jdepvar(coord,2*ii-1:2*ii),par,arcarg([1 1]),abs(jumparg(ii)));
            out2=feval(ocObj,'EventBC',jt(ii),jdepvar(coord,2*ii-1:2*ii),par,[],jumparg(ii));
            out3=feval(ocObj,'ImpulseHamiltonian',jt(ii),jdepvar(coord,2*ii-1:2*ii),par,[],jumparg(ii));
            if out1<0 || any(abs(out2)>AdmissibleTolerance)% || out3<0
                b(ii)=0;
            end
        elseif ii==numj
            out1=feval(ocObj,'InteriorImpulseBC',jt(ii),jdepvar(coord,2*ii-1:2*ii),par,arcarg([numj-1 numj-1]),abs(jumparg(ii)));
            out2=feval(ocObj,'EventBC',jt(ii),jdepvar(coord,2*ii-1:2*ii),par,[],jumparg(ii));
            out3=feval(ocObj,'ImpulseHamiltonian',jt(ii),jdepvar(coord,2*ii-1:2*ii),par,[],jumparg(ii));
            if out1>0 || any(abs(out2)>AdmissibleTolerance)% || out3<0
                b(ii)=0;
            end
        else
            out1=feval(ocObj,'InteriorImpulseBC',jt(ii),jdepvar(coord,2*ii-1:2*ii),par,arcarg([ii-1:ii]),abs(jumparg(ii)));
            out2=feval(ocObj,'EventBC',jt(ii),jdepvar(coord,2*ii-1:2*ii),par,[],jumparg(ii));
            out3=feval(ocObj,'ImpulseHamiltonian',jt(ii),jdepvar(coord,2*ii-1:2*ii),par,[],jumparg(ii));
            if any(abs([out1;out2])>AdmissibleTolerance)% || out3<0
                b(ii)=0;
            end
        end
        %val(:,ii)=[out1;out2;out3];
        val(:,ii)=[out1;out2];
    else
    end
end

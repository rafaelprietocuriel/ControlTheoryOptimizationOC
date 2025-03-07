function [b infoS newarcpos violarcarg]=testadmissibility(ocTrj,ocObj,varargin)

opt=[];
conttypeval='';
if nargin>=3
    opt=varargin{1};
end
if nargin>=4
    conttypeval=varargin{2};
end
if isempty(opt)
    opt=defaultocoptions;
end
if isempty(conttypeval)
    conttypeval=conttype(ocTrj);
end
if isempty(ocTrj)
    b=[];
    infoS=[];
    newarcpos=[];
    violarcarg=[];
    return
end
if ~isempty(conttypeval)
    switch conttypeval
        case {'extremal4fimp','extremalt4fimp','extremalp4fimp','extremalpopt4fimp'}
            modelfunc=modelspecificfunc(ocObj,'4FiniteImpulsePathContinuation');
        otherwise
            return
    end
else
    modelfunc=modelspecificfunc(ocObj,'4FiniteImpulsePathContinuation');
end
funch=modelfunc();
testadmissibility=funch{12};

AdmissibleTolerance=getocoptions(opt,'GENERAL','AdmissibleTolerance');
ZeroTimeDiffTolerance=getocoptions(opt,'GENERAL','ZeroTimeDiffTolerance');
arcarg=arcargument(ocTrj);
arcpos=arcposition(ocTrj);
switchtime=arcinterval(ocTrj);
arci=arcinterval(ocTrj);
t=time(ocObj,ocTrj,1);
t=[arci(1) t arci(end)];
modelpar=parametervalue(ocObj);
depvar=dependentvar(ocTrj);

b=0;
infoS=[];
leftarcindex=arcpos(1,:);
rightarcindex=arcpos(2,:);
counter=0;
diffarctime=diff(switchtime);
for arc=1:arcnum(ocTrj)
    constr=testadmissibility(t(leftarcindex(arc):rightarcindex(arc)),depvar(:,leftarcindex(arc):rightarcindex(arc)),modelpar,arcarg(arc),[]);
    violationmat=constr<-AdmissibleTolerance;
    if any(violationmat(:))
        counter=counter+1;
        [rows cols]=find(violationmat);
        infoS(counter).arcarg=arcarg(arc);
        infoS(counter).arcnum=arc;
        infoS(counter).rows=rows;
        infoS(counter).cols=cols;
        infoS(counter).violationmat=violationmat;
        infoS(counter).constraintvalue=constr;
        infoS(counter).minval=min(constr(:));
        b=min([b infoS(counter).minval]);
    end
    % test arctimes
    violationmat=diffarctime(arc)<-ZeroTimeDiffTolerance;
    if any(violationmat)
        counter=counter+1;
        cols=find(violationmat);
        infoS(counter).arcarg=arcarg(arc);
        infoS(counter).arcnum=arc;
        infoS(counter).rows='switchtime';
        infoS(counter).cols=cols;
        infoS(counter).violationmat=violationmat;
        infoS(counter).constraintvalue=diffarctime(arc);
        infoS(counter).minval=min(diffarctime(:));
        b=min([b infoS(counter).minval]);
    end
end

newarcpos=[];
violarcarg=[];
for ii=1:length(infoS)
    arcpos1=ocTrj.arcposition(:,infoS(ii).arcnum);
    arcpos=find(diff(infoS(ii).cols)>1).';
    newarcpos=[newarcpos arcpos1(1)-1+[infoS(ii).cols(1) infoS(ii).cols(arcpos+1).';infoS(ii).cols(arcpos).' infoS(ii).cols(end)]];
    violarcarg=[violarcarg repmat(infoS(ii).arcarg,1,length(arcpos)+1)];
end

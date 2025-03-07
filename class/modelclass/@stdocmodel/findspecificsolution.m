function [solObj,ocObjf,solidx]=findspecificsolution(ocObj,idx,varargin)
%
%

solObj=[];
ocObjf=[];
contid=[]; % c,t,i,l,m,T
value=[];
searchclass=[];
option=[];

coordinate=[]; % if value is multidimensional provide specific coordinate
if isempty(ocObj)
    return
end
contididx=find(strcmpi(varargin,'contid'));
valueidx=find(strcmpi(varargin,'value'));
coordinateidx=find(strcmpi(varargin,'coordinate'));
searchclassidx=find(strcmpi(varargin,'searchclass'));
optionidx=find(strcmpi(varargin,'option'));

if ~isempty(contididx)
    contid=varargin{contididx+1};
end
if ~isempty(valueidx)
    value=varargin{valueidx+1};
end
if ~isempty(coordinateidx)
    coordinate=varargin{coordinateidx+1};
end
if ~isempty(searchclassidx)
    searchclass=varargin{searchclassidx+1};
end
if ~isempty(optionidx)
    option=varargin{optionidx+1};
end
if isempty(coordinate)
    coordinate=1;
end
if isempty(contid)
    contid='c';
end
if isempty(value)
    value=1;
end
if isempty(searchclass)
    searchclass='exact';
end
if isempty(option)
    option=defaultocoptions;
end
if ischar(coordinate)
    coordinate=parameterindex(ocObj,coordinate);
end

contRes=contresult(ocObj);

if isempty(idx)
    idx=1:length(contRes);
end
solidx=zeros(2,length(idx));

ctr=0;
for ii=idx
    contSol=contRes{ii}.ContinuationSolution;
    contpar=zeros(length(coordinate),length(contSol));

    for jj=1:length(contSol)
        if isfield(contRes{ii},'ExtremalSolution')
            classfield=class(contRes{ii}.ExtremalSolution);
        end
        switch classfield
            case 'octrajectory'
                ocTrj=octrajectory(contSol(jj));
                switch contid
                    case 'c'
                        contpar(jj)=continuationparameter(ocTrj);
                    case 't'
                        arcint=arcinterval(ocTrj);
                        contpar(:,jj)=arcint(coordinate);
                    case 'i'
                        x=state(ocObj,ocTrj);
                        contpar(:,jj)=x(coordinate);
                    case 'l'
                        x=costate(ocObj,ocTrj);
                        contpar(:,jj)=x(coordinate);
                    case 'e'
                        x=ocTrj.y;
                        %x=state(ocObj,ocTrj,1);
                        contpar(:,jj)=x(coordinate,end);
                    case 'm'
                        par=modelparameter(ocTrj);
                        contpar(jj)=par(coordinate);
                    case 'T'
                        arcint=arcinterval(ocTrj(1));
                        contpar(jj)=arcint(end);
                end
            case 'ocasymptotic'
                ocTrj=octrajectory(contSol(jj));
                switch id
                    case 'c'
                        contpar(jj)=continuationparameter(ocTrj);
                    case 'i'
                        x=state(ocObj,ocTrj);
                        contpar(jj)=x(val(1));
                    case 'l'
                        x=costate(ocObj,ocTrj,1);
                        contpar(jj)=x(val(1));
                    case 'e'
                        x=state(ocObj,ocTrj,1);
                        contpar(jj)=x(val(1),end);
                    case 'm'
                        par=modelparameter(ocTrj);
                        contpar(jj)=par(val(1));
                end
        end
    end
    % remove identical entries
    contpar(:,contpar(1,:)==0)=[];
    [contpar,uniqueidx]=unique(contpar);
    contSol=contSol(uniqueidx);
    switch searchclass
        case 'exact'
            if length(coordinate)==1
                contidx=find(contpar==value);
            else
                contidx=find(contpar==repmat(value,1,length(contpar)));
            end
        case 'approx'
            epsilon=getocoptions(option,'GENERAL','ZeroDeviationTolerance');
            if length(coordinate)==1
                contidx=find(abs(contpar-value)<epsilon);
            else
                contidx=find(abs(contpar-repmat(value,1,length(contpar)))<repmat(epsilon,length(coordinate),length(contpar)));
            end
    end

    for kk=1:length(contidx)
        ctr=ctr+1;
        solidx(:,ctr)=[ii;contidx(kk)];
        switch classfield
            case 'octrajectory'
                solObj{ctr}=octrajectory(contSol(contidx(kk)));
                ocObjf{ctr}=changeparametervalue(ocObj,modelparameter(solObj{ctr}));

            case 'ocasymptotic'
                solObj{ctr}=octrajectory(contSol(contidx(kk)));
                ocObjf{ctr}=changeparametervalue(ocObj,modelparameter(solObj{ctr}));
                if isfield(contRes{ii},'LimitSet')
                    ocLim=contRes{ii}.LimitSet;
                else
                    eq.y=contSol(contidx(kk)).solverinfo.parameters(contSol(contidx(kk)).solverinfo.equilibriumcoord);
                    eq.arcarg=arcargument(limitset(contRes{ii}.ExtremalSolution));
                    ocLim=dynprimitive(eq);
                    ocLim.linearization=linearize(ocLim,ocObjf{ctr},'dependentvar',1);
                end
                solObj{ctr}=ocasymptotic(octrajectory(contSol(contidx(kk))),ocLim);
        end
    end
end
solidx(:,solidx(1,:)==0)=[];
function h=matcontsolutionplot(matRes,plotcoordinate)
h=[];

if isempty(matRes)
    return
end
bloodred=[175,17,28]/255;
black=[0 0 0];
blue=[0 0 1];
gray=[0.5 0.5 0.5];
lightgray=[0.75 0.75 0.75];
superlightgray=[0.95 0.95 0.95];
red=[1 0 0];
brown=[153 76 0]/255;
darkgreen=[0 0.5 0];
green=[0 1 0];
magenta=[1 0 1];
cyan=[0 1 1];
orange=[1 0.5 0];
darkviolet=[148 0 211]/255;
aquamarine=[69,139,116]/255;
rosybrown=[188,143,143]/255;
lightcoral=[240,128,128]/255;
olive=[0.5 0.5 0];
deeppink=[255,20,147]/255;
aquamarine=[127,255,212]/255;


switch matRes.ContinuationClassification
    case 'modelequilibrium'
        h=plot(matRes.ContinuationSolution.userinfo.varyparametervalue,matRes.ContinuationSolution.y(plotcoordinate,:));
    case {'modelhopf','modellimitpoint'}
        h=plot(matRes.ContinuationSolution.userinfo.varyparametervalue(1,:),matRes.ContinuationSolution.userinfo.varyparametervalue(2,:));
        for ii=1:length(matRes.ContinuationInformation)
            idx=matRes.ContinuationInformation(ii).index;
            switch matRes.ContinuationInformation(ii).label
                case 'ZH'
                  h0=plot(matRes.ContinuationSolution.userinfo.varyparametervalue(1,idx),matRes.ContinuationSolution.userinfo.varyparametervalue(2,idx));
                  set(h0,'Marker','.','MarkerSize',20,'Color',black)
            end
        end
    case 'modellimitcycle'
        actparidx=find(matRes.ContinuationInformation(end).data.parametervalues-matRes.ContinuationInformation(1).data.parametervalues);
        n=length(matRes.ContinuationSolution);
        par=zeros(1,n);
        T=par;
        for ii=1:n
            par(ii)=matRes.ContinuationSolution(ii).modelparameter(actparidx);
            T(ii)=matRes.ContinuationSolution(ii).timehorizon;
        end
        h=plot(par,T);
    case 'modellimitpointcycle'
        contSol=matRes.ContinuationSolution;
        actparidx=find(contSol(1).modelparameter-contSol(end).modelparameter);
        parval=zeros(2,length(contSol));
        for jj=1:length(contSol)
            par=contSol(jj).modelparameter;
            parval(:,jj)=par(actparidx([2 1]));
        end
        h=plot(parval(1,:),parval(2,:));
end
